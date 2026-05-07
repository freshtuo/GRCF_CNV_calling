#!/usr/bin/env python3
import argparse
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


BASE_FIELDS = [
    "chromosome",
    "segment_start",
    "segment_end",
    "segment_length",
    "log2_ratio",
    "total_cn",
    "minor_cn",
    "cnv_call",
    "loh_status",
    "event_size",
]

SEGMENT_FIELDS = BASE_FIELDS + [
    "n_genes",
    "genes",
    "max_gene_overlap_fraction",
]

GENE_FIELDS = BASE_FIELDS + [
    "gene_name",
    "gene_start",
    "gene_end",
    "gene_length_bp",
    "gene_overlap_bp",
    "gene_overlap_fraction",
    "gene_overlap_type",
]


def read_tsv(path):
    """Read a TSV file as a pandas DataFrame with string columns by default."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_genes(path):
    """Load genes_bed as a normalized BED4 table."""
    gene_frame = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
    if gene_frame.shape[1] < 3:
        raise ValueError("genes_bed must contain at least chromosome, start, and end columns")
    gene_frame = gene_frame.iloc[:, :4]
    gene_frame.columns = ["chromosome", "start", "end", "gene_name"][: gene_frame.shape[1]]
    if "gene_name" not in gene_frame.columns:
        gene_frame["gene_name"] = "."
    gene_frame["start"] = pd.to_numeric(gene_frame["start"], errors="coerce").astype("Int64")
    gene_frame["end"] = pd.to_numeric(gene_frame["end"], errors="coerce").astype("Int64")
    gene_frame = gene_frame.dropna(subset=["chromosome", "start", "end"])
    gene_frame = gene_frame.loc[gene_frame["end"] > gene_frame["start"]].copy()
    gene_frame["gene_name"] = gene_frame["gene_name"].replace("", ".")
    return gene_frame


def first_existing_column(frame, choices):
    """Return the first available column name from a list of possible names."""
    return next((name for name in choices if name in frame.columns), None)


def normalize_segments(frame, source_tool):
    """Normalize a CNVkit or FACETS segment table to shared column names."""
    column_map = {
        "chromosome": first_existing_column(frame, ["chromosome", "chrom", "chr", "Chromosome"]),
        "segment_start": first_existing_column(frame, ["start", "loc.start", "seg_start", "Start"]),
        "segment_end": first_existing_column(frame, ["end", "loc.end", "seg_end", "End"]),
        "log2_ratio": first_existing_column(frame, ["log2", "cnlr.median", "log2_ratio"]),
        "total_cn": first_existing_column(frame, ["cn", "tcn.em", "total_cn"]),
        "minor_cn": first_existing_column(frame, ["minor_cn", "lcn.em", "minor_cn.em"]),
    }
    required = ["chromosome", "segment_start", "segment_end"]
    missing = [name for name in required if column_map[name] is None]
    if missing:
        raise ValueError(f"Segment file is missing required columns: {', '.join(missing)}")

    normalized = pd.DataFrame()
    for output_col, input_col in column_map.items():
        normalized[output_col] = frame[input_col] if input_col else pd.NA

    normalized["segment_start"] = pd.to_numeric(normalized["segment_start"], errors="coerce").astype("Int64")
    normalized["segment_end"] = pd.to_numeric(normalized["segment_end"], errors="coerce").astype("Int64")
    normalized["log2_ratio"] = pd.to_numeric(normalized["log2_ratio"], errors="coerce")
    normalized["total_cn"] = pd.to_numeric(normalized["total_cn"], errors="coerce").round().astype("Int64")
    normalized["minor_cn"] = pd.to_numeric(normalized["minor_cn"], errors="coerce").round().astype("Int64")
    normalized = normalized.dropna(subset=["chromosome", "segment_start", "segment_end"]).copy()
    normalized = normalized.loc[normalized["segment_end"] > normalized["segment_start"]].copy()
    normalized.insert(0, "segment_id", [f"seg_{i + 1}" for i in range(len(normalized))])
    return normalized


def classify_segments(frame, source_tool, thresholds, focal_bp, broad_bp):
    """Add CNV class, LOH status, segment length, and event-size columns."""
    frame = frame.copy()
    frame["segment_length"] = frame["segment_end"] - frame["segment_start"]

    if source_tool == "cnvkit":
        conditions = [
            frame["log2_ratio"] > thresholds["amp"],
            frame["log2_ratio"] > thresholds["gain"],
            frame["log2_ratio"] < thresholds["deep_del"],
            frame["log2_ratio"] < thresholds["loss"],
        ]
        choices = ["amplification", "gain", "deep_deletion", "loss"]
    else:
        conditions = [
            frame["total_cn"] >= 5,
            frame["total_cn"] >= 3,
            frame["total_cn"] <= 0,
            frame["total_cn"] <= 1,
        ]
        choices = ["amplification", "gain", "homozygous_deletion", "loss"]

    frame["cnv_call"] = pd.NA
    for condition, choice in zip(conditions, choices):
        frame.loc[condition & frame["cnv_call"].isna(), "cnv_call"] = choice
    frame["cnv_call"] = frame["cnv_call"].fillna("neutral")
    missing_call = frame["log2_ratio"].isna() if source_tool == "cnvkit" else frame["total_cn"].isna()
    frame.loc[missing_call, "cnv_call"] = "unknown"

    frame["loh_status"] = "retained"
    frame.loc[frame["minor_cn"] == 0, "loh_status"] = "LOH"
    frame["event_size"] = "intermediate"
    frame.loc[frame["segment_length"] < focal_bp, "event_size"] = "focal"
    frame.loc[frame["segment_length"] >= broad_bp, "event_size"] = "broad"
    return frame


def write_bedtools_inputs(segments, genes, segment_bed, gene_bed):
    """Write BED files used by bedtools intersect."""
    segments[["chromosome", "segment_start", "segment_end", "segment_id"]].to_csv(
        segment_bed, sep="\t", header=False, index=False
    )
    genes[["chromosome", "start", "end", "gene_name"]].to_csv(
        gene_bed, sep="\t", header=False, index=False
    )


def run_bedtools_intersect(segments, genes):
    """Return one row per segment-gene overlap using bedtools intersect -wo."""
    if segments.empty or genes.empty:
        return pd.DataFrame(
            columns=[
                "segment_chromosome",
                "segment_start",
                "segment_end",
                "segment_id",
                "gene_chromosome",
                "gene_start",
                "gene_end",
                "gene_name",
                "gene_overlap_bp",
            ]
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        segment_bed = Path(tmpdir) / "segments.bed"
        gene_bed = Path(tmpdir) / "genes.bed"
        write_bedtools_inputs(segments, genes, segment_bed, gene_bed)
        command = [
            "bedtools",
            "intersect",
            "-a",
            str(segment_bed),
            "-b",
            str(gene_bed),
            "-wo",
        ]
        result = subprocess.run(command, check=True, capture_output=True, text=True)

    if not result.stdout.strip():
        return pd.DataFrame(
            columns=[
                "segment_chromosome",
                "segment_start",
                "segment_end",
                "segment_id",
                "gene_chromosome",
                "gene_start",
                "gene_end",
                "gene_name",
                "gene_overlap_bp",
            ]
        )

    rows = [line.split("\t") for line in result.stdout.rstrip("\n").splitlines()]
    overlaps = pd.DataFrame(
        rows,
        columns=[
            "segment_chromosome",
            "segment_start",
            "segment_end",
            "segment_id",
            "gene_chromosome",
            "gene_start",
            "gene_end",
            "gene_name",
            "gene_overlap_bp",
        ],
    )
    for column in ["segment_start", "segment_end", "gene_start", "gene_end", "gene_overlap_bp"]:
        overlaps[column] = pd.to_numeric(overlaps[column], errors="coerce").astype("Int64")
    return overlaps


def annotate_gene_overlaps(segments, genes, min_overlap_bp, min_overlap_fraction):
    """Intersect segments and genes, then calculate overlap metrics."""
    overlaps = run_bedtools_intersect(segments, genes)
    if overlaps.empty:
        return overlaps

    overlaps["gene_length_bp"] = overlaps["gene_end"] - overlaps["gene_start"]
    overlaps["gene_overlap_fraction"] = (
        overlaps["gene_overlap_bp"].astype(float) / overlaps["gene_length_bp"].astype(float)
    )
    overlaps["gene_overlap_type"] = "partial"
    overlaps.loc[
        (overlaps["gene_overlap_bp"] == overlaps["gene_length_bp"])
        & (overlaps["gene_length_bp"] > 0),
        "gene_overlap_type",
    ] = "complete"
    overlaps = overlaps.loc[
        (overlaps["gene_overlap_bp"] >= min_overlap_bp)
        & (overlaps["gene_overlap_fraction"] >= min_overlap_fraction)
    ].copy()
    return overlaps


def build_segment_output(segments, overlaps):
    """Write one row per segment with collapsed gene annotations."""
    output = segments.copy()
    if overlaps.empty:
        output["n_genes"] = 0
        output["genes"] = ""
        output["max_gene_overlap_fraction"] = ""
    else:
        summary = (
            overlaps.sort_values(["segment_id", "gene_name"])
            .groupby("segment_id")
            .agg(
                n_genes=("gene_name", "nunique"),
                genes=("gene_name", lambda values: ",".join(sorted(set(v for v in values if v and v != ".")))),
                max_gene_overlap_fraction=("gene_overlap_fraction", "max"),
            )
            .reset_index()
        )
        output = output.merge(summary, on="segment_id", how="left")
        output["n_genes"] = output["n_genes"].fillna(0).astype(int)
        output["genes"] = output["genes"].fillna("")
        output["max_gene_overlap_fraction"] = output["max_gene_overlap_fraction"].fillna("")
    return output[SEGMENT_FIELDS]


def build_gene_output(segments, overlaps):
    """Write one compact row per affected gene."""
    if overlaps.empty:
        return pd.DataFrame(columns=GENE_FIELDS)

    overlap_columns = overlaps.drop(columns=["segment_chromosome", "segment_start", "segment_end"])
    gene_rows = overlap_columns.merge(segments, on="segment_id", how="left")
    gene_rows = gene_rows.loc[~gene_rows["cnv_call"].isin(["neutral", "unknown"])].copy()
    if gene_rows.empty:
        return pd.DataFrame(columns=GENE_FIELDS)
    gene_rows = gene_rows.sort_values(
        ["gene_name", "gene_overlap_fraction", "gene_overlap_bp", "segment_length"],
        ascending=[True, False, False, False],
    )
    gene_rows = gene_rows.drop_duplicates(subset=["gene_name"], keep="first")
    return gene_rows[
        BASE_FIELDS
        + [
            "gene_name",
            "gene_start",
            "gene_end",
            "gene_length_bp",
            "gene_overlap_bp",
            "gene_overlap_fraction",
            "gene_overlap_type",
        ]
    ]


def main():
    """Annotate one comparison's CNV segments and write segment/gene TSV outputs."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-tool", choices=["cnvkit", "facets"], required=True)
    parser.add_argument("--segments", required=True)
    parser.add_argument("--genes-bed", required=True)
    parser.add_argument("--thresholds", required=True, help="gain,loss,amp,deep_del,focal_mb,broad_mb")
    parser.add_argument("--min-gene-overlap-bp", type=int, default=1)
    parser.add_argument("--min-gene-overlap-fraction", type=float, default=0.0)
    parser.add_argument("--out-segments", required=True)
    parser.add_argument("--out-genes", required=True)
    args = parser.parse_args()

    # Thresholds arrive as one comma-separated Snakemake parameter so the rule
    # can stay simple while still keeping classification settings in config.yaml.
    gain, loss, amp, deep_del, focal_mb, broad_mb = [float(x) for x in args.thresholds.split(",")]
    thresholds = {"gain": gain, "loss": loss, "amp": amp, "deep_del": deep_del}
    focal_bp = int(focal_mb * 1_000_000)
    broad_bp = int(broad_mb * 1_000_000)

    genes = load_genes(args.genes_bed)
    segment_rows = classify_segments(
        normalize_segments(read_tsv(args.segments), args.source_tool),
        args.source_tool,
        thresholds,
        focal_bp,
        broad_bp,
    )
    overlaps = annotate_gene_overlaps(
        segment_rows,
        genes,
        args.min_gene_overlap_bp,
        args.min_gene_overlap_fraction,
    )

    Path(args.out_segments).parent.mkdir(parents=True, exist_ok=True)
    build_segment_output(segment_rows, overlaps).to_csv(args.out_segments, sep="\t", index=False, na_rep="")
    build_gene_output(segment_rows, overlaps).to_csv(args.out_genes, sep="\t", index=False, na_rep="")


if __name__ == "__main__":
    main()
