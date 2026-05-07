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


def read_optional_tsv(path):
    """Read a TSV file if it exists."""
    if not Path(path).exists():
        return pd.DataFrame()
    return read_tsv(path)


def caller_paths(results_dir, comparison_id, caller):
    """Return detailed annotation paths for one caller."""
    base = Path(results_dir) / caller / comparison_id / "annotation"
    return {
        "segments": base / "annotated_segments.tsv",
        "genes": base / "annotated_genes.tsv",
    }


def detail_links(results_dir, comparison, caller):
    """Render links to detailed per-caller outputs."""
    cid = comparison["comparison_id"]
    paths = caller_paths(results_dir, cid, caller)
    if not paths["segments"].exists() and not paths["genes"].exists():
        return ""
    rel_segments = f"../../{caller}/{cid}/annotation/annotated_segments.tsv"
    rel_genes = f"../../{caller}/{cid}/annotation/annotated_genes.tsv"
    return (
        f"<li>{escape(caller)}: "
        f"<a href=\"{escape(rel_segments)}\">annotated segments</a>, "
        f"<a href=\"{escape(rel_genes)}\">annotated genes</a></li>"
    )


def top_genes(results_dir, comparison):
    """Collect top affected-gene rows from each enabled caller."""
    frames = []
    cid = comparison["comparison_id"]
    for caller in ("cnvkit", "facets"):
        if not yes(comparison.get(f"run_{caller}", "")):
            continue
        genes = read_optional_tsv(caller_paths(results_dir, cid, caller)["genes"])
        if genes.empty:
            continue
        genes.insert(0, "caller", caller)
        frames.append(genes)
    if not frames:
        return pd.DataFrame()
    frame = pd.concat(frames, ignore_index=True, sort=False)
    preferred = [
        "caller",
        "gene_name",
        "cnv_call",
        "loh_status",
        "event_size",
        "chromosome",
        "segment_start",
        "segment_end",
        "gene_overlap_fraction",
    ]
    columns = [column for column in preferred if column in frame.columns]
    return frame[columns]


def write_comparison_report(report_path, comparison, qc, cnv_summary, purity_ploidy, results_dir):
    """Write one detailed HTML report for a comparison."""
    cid = comparison["comparison_id"]
    qc_frame = qc.loc[qc["comparison_id"] == cid].copy() if not qc.empty else pd.DataFrame()
    summary_frame = (
        cnv_summary.loc[cnv_summary["comparison_id"] == cid].copy()
        if not cnv_summary.empty
        else pd.DataFrame()
    )
    pp_frame = (
        purity_ploidy.loc[purity_ploidy["comparison_id"] == cid].copy()
        if not purity_ploidy.empty
        else pd.DataFrame()
    )
    genes_frame = top_genes(results_dir, comparison)
    control = comparison.get("control_id", "").strip() or "none"
    caution = ""
    if comparison.get("comparison_type", "") == "tumor_only" or not comparison.get("control_id", "").strip():
        caution = "<p><strong>Note:</strong> Tumor-only CNV calls should be interpreted cautiously because germline CNVs cannot be subtracted.</p>"

    links = []
    for caller in ("cnvkit", "facets"):
        if yes(comparison.get(f"run_{caller}", "")):
            link = detail_links(results_dir, comparison, caller)
            if link:
                links.append(link)
    detail_section = "<ul>" + "".join(links) + "</ul>" if links else "<p>No detailed annotation files available.</p>"

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
    <p><strong>FACETS:</strong> {"run" if yes(comparison.get("run_facets", "")) else "not run"}</p>
  </div>
  {caution}
  <h2>QC</h2>
  {html_table(qc_frame)}
  <h2>Purity And Ploidy</h2>
  {html_table(pp_frame)}
  <h2>CNV Summary</h2>
  {html_table(summary_frame)}
  <h2>Top Affected Genes</h2>
  {html_table(genes_frame)}
  <h2>Detailed Files</h2>
  {detail_section}
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
    parser.add_argument("--cnv-summary", required=True)
    parser.add_argument("--purity-ploidy", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--reports-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    comparisons = read_tsv(args.comparisons)
    qc = read_tsv(args.qc)
    cnv_summary = read_tsv(args.cnv_summary)
    purity_ploidy = read_tsv(args.purity_ploidy)
    reports_dir = Path(args.reports_dir)
    reports_dir.mkdir(parents=True, exist_ok=True)

    links = []
    for comparison in comparisons.to_dict(orient="records"):
        cid = comparison["comparison_id"]
        report_name = f"{cid}.report.html"
        write_comparison_report(
            reports_dir / report_name,
            comparison,
            qc,
            cnv_summary,
            purity_ploidy,
            args.results_dir,
        )
        rows = cnv_summary.loc[cnv_summary["comparison_id"] == cid] if not cnv_summary.empty else pd.DataFrame()
        callers = ",".join(rows["caller"].tolist()) if not rows.empty and "caller" in rows.columns else ""
        affected = rows["n_affected_genes"].astype(int).sum() if not rows.empty and "n_affected_genes" in rows.columns else 0
        links.append(
            "<tr>"
            f"<td><a href=\"reports/{escape(report_name)}\">{escape(cid)}</a></td>"
            f"<td>{escape(comparison['patient_id'])}</td>"
            f"<td>{escape(comparison['comparison_type'])}</td>"
            f"<td>{escape(callers)}</td>"
            f"<td>{affected}</td>"
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
    <thead><tr><th>Comparison</th><th>Patient</th><th>Type</th><th>Callers</th><th>Affected Genes</th></tr></thead>
    <tbody>{''.join(links)}</tbody>
  </table>
  <h2>Merged Tables</h2>
  <ul>
    <li><a href="all_comparisons.qc.tsv">all_comparisons.qc.tsv</a></li>
    <li><a href="all_comparisons.cnv_summary.tsv">all_comparisons.cnv_summary.tsv</a></li>
    <li><a href="all_comparisons.purity_ploidy.tsv">all_comparisons.purity_ploidy.tsv</a></li>
  </ul>
</body>
</html>
"""
    with open(args.output, "w") as handle:
        handle.write(index_html)


if __name__ == "__main__":
    main()
