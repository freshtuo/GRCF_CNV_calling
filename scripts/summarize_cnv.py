#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def read_tsv(path):
    """Read a TSV file as a pandas DataFrame with blanks preserved."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def parse_flagstat(path):
    """Extract a few high-level read-count metrics from samtools flagstat output."""
    summary = {"total_reads": "", "mapped_reads": "", "properly_paired_reads": ""}
    if not Path(path).exists():
        return summary
    with open(path) as handle:
        for line in handle:
            count = line.split(" + ")[0]
            if " in total " in line:
                summary["total_reads"] = count
            elif " mapped (" in line and "primary mapped" not in line:
                summary["mapped_reads"] = count
            elif " properly paired (" in line:
                summary["properly_paired_reads"] = count
    return summary


def comparison_qc_rows(samples, comparisons, results):
    """Build one QC summary row per comparison/sample role."""
    sample_by_id = samples.set_index("sample_id", drop=False).to_dict(orient="index")
    rows = []
    for comparison in comparisons.to_dict(orient="records"):
        cid = comparison["comparison_id"]
        sample_roles = [("case", comparison["case_id"])]
        control_id = comparison.get("control_id", "").strip()
        if control_id:
            sample_roles.append(("control", control_id))

        for role, sample_id in sample_roles:
            sample = sample_by_id.get(sample_id, {})
            flagstat = parse_flagstat(f"{results}/qc/{cid}/flagstat.{sample_id}.txt")
            quickcheck = Path(f"{results}/qc/{cid}/quickcheck.{sample_id}.txt")
            status = "pass" if quickcheck.exists() and quickcheck.stat().st_size == 0 else "review"
            rows.append(
                {
                    "comparison_id": cid,
                    "patient_id": comparison["patient_id"],
                    "sample_role": role,
                    "sample_id": sample_id,
                    "species": sample.get("species", ""),
                    "sample_type": sample.get("sample_type", ""),
                    "bam": sample.get("bam", ""),
                    **flagstat,
                    "quickcheck_status": status,
                }
            )
    return rows


def add_context(frame, comparison, caller):
    """Add comparison and caller metadata to an annotated CNV table."""
    frame = frame.copy()
    frame.insert(0, "caller", caller)
    frame.insert(0, "comparison_type", comparison["comparison_type"])
    frame.insert(0, "species", comparison["species"])
    frame.insert(0, "patient_id", comparison["patient_id"])
    frame.insert(0, "comparison_id", comparison["comparison_id"])
    return frame


def read_optional_tsv(path):
    """Read a TSV file if it exists, otherwise return an empty table."""
    if not Path(path).exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def merged_annotation_rows(comparisons, results, table_name):
    """Merge caller annotation tables across all enabled comparisons."""
    frames = []
    for comparison in comparisons.to_dict(orient="records"):
        cid = comparison["comparison_id"]
        if str(comparison.get("run_cnvkit", "")).strip().lower() in {"yes", "true", "1", "y"}:
            path = f"{results}/cnvkit/{cid}/annotation/annotated_{table_name}.tsv"
            frame = read_optional_tsv(path)
            if not frame.empty:
                frames.append(add_context(frame, comparison, "cnvkit"))
        if str(comparison.get("run_facets", "")).strip().lower() in {"yes", "true", "1", "y"}:
            path = f"{results}/facets/{cid}/annotation/annotated_{table_name}.tsv"
            frame = read_optional_tsv(path)
            if not frame.empty:
                frames.append(add_context(frame, comparison, "facets"))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def main():
    """Build final project summary tables from per-sample and per-comparison outputs."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--comparisons", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--out-qc", required=True)
    parser.add_argument("--out-segments", required=True)
    parser.add_argument("--out-genes", required=True)
    args = parser.parse_args()

    samples = read_tsv(args.samples)
    comparisons = read_tsv(args.comparisons)

    pd.DataFrame(comparison_qc_rows(samples, comparisons, args.results)).to_csv(
        args.out_qc, sep="\t", index=False
    )
    merged_annotation_rows(comparisons, args.results, "segments").to_csv(
        args.out_segments, sep="\t", index=False
    )
    merged_annotation_rows(comparisons, args.results, "genes").to_csv(
        args.out_genes, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
