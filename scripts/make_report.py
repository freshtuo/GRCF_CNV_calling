#!/usr/bin/env python3
import argparse
from html import escape

import pandas as pd


def read_tsv(path):
    """Read a TSV file, returning an empty DataFrame for empty summary files."""
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def count_rows(path):
    return len(read_tsv(path))


def count_rows_many(paths):
    return sum(count_rows(path) for path in paths)


def call_counts(paths):
    counts = []
    for path in paths:
        segments = read_tsv(path)
        if "cnv_call" in segments.columns:
            counts.append(segments["cnv_call"].fillna("unknown").value_counts())
    if not counts:
        return pd.Series(dtype=int)
    return pd.concat(counts, axis=1).fillna(0).sum(axis=1).astype(int).sort_index()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--cnvkit-segments", nargs="*", default=[])
    parser.add_argument("--cnvkit-genes", nargs="*", default=[])
    parser.add_argument("--facets-segments", nargs="*", default=[])
    parser.add_argument("--facets-genes", nargs="*", default=[])
    parser.add_argument("--qc", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    counts = call_counts(args.cnvkit_segments + args.facets_segments)
    rows = "\n".join(f"<tr><td>{escape(call)}</td><td>{count}</td></tr>" for call, count in counts.items())
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{escape(args.project)} CNV report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
    table {{ border-collapse: collapse; margin-top: 1rem; }}
    th, td {{ border: 1px solid #ccc; padding: 0.4rem 0.7rem; text-align: left; }}
  </style>
</head>
<body>
  <h1>{escape(args.project)} CNV report</h1>
  <p>CNVkit segment rows: {count_rows_many(args.cnvkit_segments)}</p>
  <p>CNVkit affected gene rows: {count_rows_many(args.cnvkit_genes)}</p>
  <p>FACETS segment rows: {count_rows_many(args.facets_segments)}</p>
  <p>FACETS affected gene rows: {count_rows_many(args.facets_genes)}</p>
  <p>QC sample rows: {count_rows(args.qc)}</p>
  <h2>CNV calls</h2>
  <table>
    <thead><tr><th>Call</th><th>Rows</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <p>Tumor-only CNV calls should be interpreted cautiously because germline CNVs cannot be removed without a matched control.</p>
</body>
</html>
"""
    with open(args.output, "w") as handle:
        handle.write(html)


if __name__ == "__main__":
    main()
