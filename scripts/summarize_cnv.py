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


def main():
    """Build final project summary tables from per-sample and per-comparison outputs."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--out-qc", required=True)
    args = parser.parse_args()

    samples = read_tsv(args.samples)

    # QC summary is sample-oriented because BAM QC runs once per physical sample.
    qc_rows = []
    for sample in samples.to_dict(orient="records"):
        sid = sample["sample_id"]
        flagstat = parse_flagstat(f"{args.results}/qc/samples/{sid}/flagstat.txt")
        quickcheck = Path(f"{args.results}/qc/samples/{sid}/quickcheck.txt")
        # samtools quickcheck is quiet on success; non-empty output needs review.
        status = "pass" if quickcheck.exists() and quickcheck.stat().st_size == 0 else "review"
        qc_rows.append(
            {
                "sample_id": sid,
                "patient_id": sample["patient_id"],
                "species": sample["species"],
                "sample_type": sample["sample_type"],
                **flagstat,
                "quickcheck_status": status,
            }
        )
    pd.DataFrame(qc_rows).to_csv(args.out_qc, sep="\t", index=False)


if __name__ == "__main__":
    main()
