# Merge comparison-scoped QC and caller-level CNV counts into project summaries.
# Detailed annotations stay in per-comparison, per-caller output folders.
rule summarize_results:
    input:
        qc_flagstat=expand(
            f"{RESULTS}/qc/{{comparison_id}}/flagstat.{{sample_id}}.txt",
            zip,
            comparison_id=QC_COMPARISON_IDS,
            sample_id=QC_SAMPLE_IDS,
        ),
        qc_quickcheck=expand(
            f"{RESULTS}/qc/{{comparison_id}}/quickcheck.{{sample_id}}.txt",
            zip,
            comparison_id=QC_COMPARISON_IDS,
            sample_id=QC_SAMPLE_IDS,
        ),
        cnvkit_segments=expand(CNVKIT_ANNOTATED_SEGMENTS, comparison_id=CNVKIT_COMPARISONS),
        cnvkit_genes=expand(CNVKIT_ANNOTATED_GENES, comparison_id=CNVKIT_COMPARISONS),
        facets_segments=expand(FACETS_ANNOTATED_SEGMENTS, comparison_id=FACETS_COMPARISONS),
        facets_genes=expand(FACETS_ANNOTATED_GENES, comparison_id=FACETS_COMPARISONS)
    output:
        qc=f"{RESULTS}/summary/all_comparisons.qc.tsv",
        cnv_summary=f"{RESULTS}/summary/all_comparisons.cnv_summary.tsv",
        purity_ploidy=f"{RESULTS}/summary/all_comparisons.purity_ploidy.tsv"
    log:
        f"{RESULTS}/summary/summarize.log"
    params:
        samples=SAMPLES_TSV,
        comparisons=COMPARISONS_TSV,
        results=lambda wildcards, output: str(Path(output.qc).parents[1])
    conda:
        "../../envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {log}); "
        "python scripts/summarize_cnv.py "
        "--samples {params.samples} "
        "--comparisons {params.comparisons} "
        "--results {params.results} "
        "--out-qc {output.qc} "
        "--out-cnv-summary {output.cnv_summary} "
        "--out-purity-ploidy {output.purity_ploidy} > {log} 2>&1"


# Create a lightweight HTML report from separate caller annotations and QC.
rule make_report:
    input:
        qc=f"{RESULTS}/summary/all_comparisons.qc.tsv",
        cnv_summary=f"{RESULTS}/summary/all_comparisons.cnv_summary.tsv",
        purity_ploidy=f"{RESULTS}/summary/all_comparisons.purity_ploidy.tsv"
    output:
        index=f"{RESULTS}/summary/report.html",
        comparison_reports=expand(f"{RESULTS}/summary/reports/{{comparison_id}}.report.html", comparison_id=COMPARISON_IDS)
    log:
        f"{RESULTS}/summary/report.log"
    params:
        project=PROJECT,
        comparisons=COMPARISONS_TSV,
        results_dir=lambda wildcards, output: str(Path(output.index).parents[1]),
        reports_dir=lambda wildcards, output: str(Path(output.index).parent / "reports")
    conda:
        "../../envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {log}); "
        "python scripts/make_report.py "
        "--project {params.project} "
        "--comparisons {params.comparisons} "
        "--qc {input.qc} "
        "--cnv-summary {input.cnv_summary} "
        "--purity-ploidy {input.purity_ploidy} "
        "--results-dir {params.results_dir} "
        "--reports-dir {params.reports_dir} "
        "--output {output.index} > {log} 2>&1"
