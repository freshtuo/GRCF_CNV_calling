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
                    "bam_name": Path(sample.get("bam", "")).name,
                    **flagstat,
                    "quickcheck_status": status,
                }
            )
    return rows


def read_optional_tsv(path):
    """Read a TSV file if it exists, otherwise return an empty table."""
    if not Path(path).exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def count_value(frame, column, value):
    """Count rows where a column equals a value."""
    if frame.empty or column not in frame.columns:
        return 0
    return int((frame[column] == value).sum())


def affected_gene_count(frame):
    """Count unique genes with non-neutral CNV calls or LOH."""
    if frame.empty or "gene_name" not in frame.columns:
        return 0
    genes = frame["gene_name"].fillna("").astype(str)
    affected = pd.Series(False, index=frame.index)
    if "cnv_call" in frame.columns:
        affected |= ~frame["cnv_call"].fillna("").isin(["", "neutral", "unknown"])
    if "loh_status" in frame.columns:
        affected |= frame["loh_status"].fillna("") == "LOH"
    genes = genes.loc[affected & (genes != "") & (genes != ".")]
    return int(genes.nunique())


def purity_ploidy_path(results, comparison_id):
    """Return the FACETS purity/ploidy output path for a comparison."""
    return f"{results}/facets/{comparison_id}/purity_ploidy/facets_purity_ploidy.tsv"


def read_purity_ploidy(results, comparison_id):
    """Read FACETS purity/ploidy values if available."""
    frame = read_optional_tsv(purity_ploidy_path(results, comparison_id))
    if frame.empty:
        return {"purity": "", "ploidy": "", "purity_status": "not_available"}
    row = frame.iloc[0].to_dict()
    purity = row.get("purity", "")
    ploidy = row.get("ploidy", "")
    status = "estimated" if str(purity).strip() else "not_estimated"
    return {"purity": purity, "ploidy": ploidy, "purity_status": status}


def cnv_summary_rows(comparisons, results):
    """Build one CNV overview row per comparison and caller."""
    rows = []
    for comparison in comparisons.to_dict(orient="records"):
        cid = comparison["comparison_id"]
        if str(comparison.get("run_cnvkit", "")).strip().lower() in {"yes", "true", "1", "y"}:
            segments = read_optional_tsv(f"{results}/cnvkit/{cid}/annotation/annotated_segments.tsv")
            genes = read_optional_tsv(f"{results}/cnvkit/{cid}/annotation/annotated_genes.tsv")
            rows.append(summary_row(comparison, "cnvkit", segments, genes, "", ""))
        if str(comparison.get("run_facets", "")).strip().lower() in {"yes", "true", "1", "y"}:
            segments = read_optional_tsv(f"{results}/facets/{cid}/annotation/annotated_segments.tsv")
            genes = read_optional_tsv(f"{results}/facets/{cid}/annotation/annotated_genes.tsv")
            pp = read_purity_ploidy(results, cid)
            rows.append(
                summary_row(
                    comparison,
                    "facets",
                    segments,
                    genes,
                    pp["purity"],
                    pp["ploidy"],
                    pp["purity_status"],
                )
            )
    return rows


def summary_row(comparison, caller, segments, genes, purity, ploidy, purity_status=""):
    """Summarize CNV calls for one comparison/caller."""
    return {
        "comparison_id": comparison["comparison_id"],
        "patient_id": comparison["patient_id"],
        "species": comparison["species"],
        "comparison_type": comparison["comparison_type"],
        "caller": caller,
        "n_segments": len(segments),
        "n_gain": count_value(segments, "cnv_call", "gain"),
        "n_loss": count_value(segments, "cnv_call", "loss"),
        "n_amplification": count_value(segments, "cnv_call", "amplification"),
        "n_deep_deletion": count_value(segments, "cnv_call", "deep_deletion"),
        "n_homozygous_deletion": count_value(segments, "cnv_call", "homozygous_deletion"),
        "n_loh_segments": count_value(segments, "loh_status", "LOH"),
        "n_focal_segments": count_value(segments, "event_size", "focal"),
        "n_broad_segments": count_value(segments, "event_size", "broad"),
        "n_affected_genes": affected_gene_count(genes),
        "purity": purity,
        "ploidy": ploidy,
        "purity_status": purity_status,
    }


def purity_ploidy_rows(comparisons, results):
    """Build one purity/ploidy row per comparison, blank when unavailable."""
    rows = []
    for comparison in comparisons.to_dict(orient="records"):
        pp = read_purity_ploidy(results, comparison["comparison_id"])
        rows.append(
            {
                "comparison_id": comparison["comparison_id"],
                "patient_id": comparison["patient_id"],
                "comparison_type": comparison["comparison_type"],
                "purity": pp["purity"],
                "ploidy": pp["ploidy"],
                "purity_status": pp["purity_status"],
            }
        )
    return rows


def main():
    """Build final project summary tables from per-sample and per-comparison outputs."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--comparisons", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--out-qc", required=True)
    parser.add_argument("--out-cnv-summary", required=True)
    parser.add_argument("--out-purity-ploidy", required=True)
    args = parser.parse_args()

    samples = read_tsv(args.samples)
    comparisons = read_tsv(args.comparisons)

    pd.DataFrame(comparison_qc_rows(samples, comparisons, args.results)).to_csv(
        args.out_qc, sep="\t", index=False
    )
    pd.DataFrame(cnv_summary_rows(comparisons, args.results)).to_csv(
        args.out_cnv_summary, sep="\t", index=False
    )
    pd.DataFrame(purity_ploidy_rows(comparisons, args.results)).to_csv(
        args.out_purity_ploidy, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
