#!/usr/bin/env python3
import argparse
from html import escape
from pathlib import Path

import pandas as pd


def read_tsv(path):
    """Read a TSV file, returning an empty DataFrame for empty files."""
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def yes(value):
    """Normalize yes/no style metadata fields to a boolean decision."""
    return str(value).strip().lower() in {"yes", "true", "1", "y"}


def html_table(frame, max_rows=50):
    """Render a compact HTML table for a dataframe."""
    if frame.empty:
        return "<p>No rows.</p>"
    display = frame.head(max_rows)
    headers = "".join(f"<th>{escape(str(col))}</th>" for col in display.columns)
    body_rows = []
    for row in display.to_dict(orient="records"):
        cells = "".join(f"<td>{escape(str(value))}</td>" for value in row.values())
        body_rows.append(f"<tr>{cells}</tr>")
    note = ""
    if len(frame) > max_rows:
        note = f"<p>Showing first {max_rows} of {len(frame)} rows.</p>"
    return f"{note}<table><thead><tr>{headers}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def call_counts(frame):
    """Return CNV call counts from a segment table."""
    if frame.empty or "cnv_call" not in frame.columns:
        return pd.DataFrame(columns=["cnv_call", "rows"])
    counts = frame["cnv_call"].replace("", "unknown").value_counts().sort_index()
    return counts.rename_axis("cnv_call").reset_index(name="rows")


def write_comparison_report(report_path, comparison, qc, segments, genes):
    """Write one detailed HTML report for a comparison."""
    cid = comparison["comparison_id"]
    qc_frame = qc.loc[qc["comparison_id"] == cid].copy() if not qc.empty else pd.DataFrame()
    segment_frame = segments.loc[segments["comparison_id"] == cid].copy() if not segments.empty else pd.DataFrame()
    gene_frame = genes.loc[genes["comparison_id"] == cid].copy() if not genes.empty else pd.DataFrame()
    counts = call_counts(segment_frame)
    facets_status = "run" if yes(comparison.get("run_facets", "")) else "not run"
    control = comparison.get("control_id", "").strip() or "none"
    caution = ""
    if comparison.get("comparison_type", "") == "tumor_only" or not comparison.get("control_id", "").strip():
        caution = "<p><strong>Note:</strong> Tumor-only CNV calls should be interpreted cautiously because germline CNVs cannot be subtracted.</p>"

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{escape(cid)} CNV report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
    h1, h2 {{ margin-top: 1.4rem; }}
    table {{ border-collapse: collapse; margin-top: 0.8rem; font-size: 0.9rem; }}
    th, td {{ border: 1px solid #ccc; padding: 0.35rem 0.55rem; text-align: left; vertical-align: top; }}
    th {{ background: #f5f5f5; }}
    .meta p {{ margin: 0.2rem 0; }}
  </style>
</head>
<body>
  <p><a href="../report.html">Project index</a></p>
  <h1>{escape(cid)}</h1>
  <div class="meta">
    <p><strong>Patient:</strong> {escape(comparison["patient_id"])}</p>
    <p><strong>Species:</strong> {escape(comparison["species"])}</p>
    <p><strong>Type:</strong> {escape(comparison["comparison_type"])}</p>
    <p><strong>Case:</strong> {escape(comparison["case_id"])}</p>
    <p><strong>Control:</strong> {escape(control)}</p>
    <p><strong>CNVkit:</strong> {"run" if yes(comparison.get("run_cnvkit", "")) else "not run"}</p>
    <p><strong>FACETS:</strong> {escape(facets_status)}</p>
  </div>
  {caution}
  <h2>QC</h2>
  {html_table(qc_frame)}
  <h2>CNV Call Counts</h2>
  {html_table(counts)}
  <h2>Affected Genes</h2>
  {html_table(gene_frame)}
  <h2>Segments</h2>
  {html_table(segment_frame)}
</body>
</html>
"""
    with open(report_path, "w") as handle:
        handle.write(html)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--comparisons", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--segments", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--reports-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    comparisons = read_tsv(args.comparisons)
    qc = read_tsv(args.qc)
    segments = read_tsv(args.segments)
    genes = read_tsv(args.genes)
    reports_dir = Path(args.reports_dir)
    reports_dir.mkdir(parents=True, exist_ok=True)

    links = []
    for comparison in comparisons.to_dict(orient="records"):
        cid = comparison["comparison_id"]
        report_name = f"{cid}.report.html"
        write_comparison_report(reports_dir / report_name, comparison, qc, segments, genes)
        segment_count = len(segments.loc[segments["comparison_id"] == cid]) if not segments.empty else 0
        gene_count = len(genes.loc[genes["comparison_id"] == cid]) if not genes.empty else 0
        links.append(
            "<tr>"
            f"<td><a href=\"reports/{escape(report_name)}\">{escape(cid)}</a></td>"
            f"<td>{escape(comparison['patient_id'])}</td>"
            f"<td>{escape(comparison['comparison_type'])}</td>"
            f"<td>{segment_count}</td>"
            f"<td>{gene_count}</td>"
            "</tr>"
        )

    index_html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{escape(args.project)} CNV report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
    table {{ border-collapse: collapse; margin-top: 1rem; }}
    th, td {{ border: 1px solid #ccc; padding: 0.4rem 0.7rem; text-align: left; }}
    th {{ background: #f5f5f5; }}
  </style>
</head>
<body>
  <h1>{escape(args.project)} CNV report</h1>
  <h2>Comparisons</h2>
  <table>
    <thead><tr><th>Comparison</th><th>Patient</th><th>Type</th><th>Segment Rows</th><th>Gene Rows</th></tr></thead>
    <tbody>{''.join(links)}</tbody>
  </table>
  <h2>Merged Tables</h2>
  <ul>
    <li><a href="all_comparisons.qc.tsv">all_comparisons.qc.tsv</a></li>
    <li><a href="all_comparisons.segments.tsv">all_comparisons.segments.tsv</a></li>
    <li><a href="all_comparisons.genes.tsv">all_comparisons.genes.tsv</a></li>
  </ul>
</body>
</html>
"""
    with open(args.output, "w") as handle:
        handle.write(index_html)


if __name__ == "__main__":
    main()
