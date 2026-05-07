# Merge per-sample QC into one project table. CNV caller annotations stay in
# their per-comparison, per-caller files to keep interpretation contexts clear.
rule summarize_results:
    input:
        qc=expand(f"{RESULTS}/qc/samples/{{sample_id}}/flagstat.txt", sample_id=SAMPLES.keys())
    output:
        qc=f"{RESULTS}/summary/qc.tsv"
    log:
        f"{RESULTS}/summary/summarize.log"
    params:
        samples=SAMPLES_TSV,
        results=lambda wildcards, output: str(Path(output.qc).parents[1])
    conda:
        "envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {log}); "
        "python scripts/summarize_cnv.py "
        "--samples {params.samples} "
        "--results {params.results} "
        "--out-qc {output.qc} > {log} 2>&1"


# Create a lightweight HTML report from separate caller annotations and QC.
rule make_report:
    input:
        cnvkit_segments=expand(f"{RESULTS}/annotation/{{comparison_id}}/annotated_segments.tsv", comparison_id=CNVKIT_COMPARISONS),
        cnvkit_genes=expand(f"{RESULTS}/annotation/{{comparison_id}}/annotated_genes.tsv", comparison_id=CNVKIT_COMPARISONS),
        facets_segments=expand(f"{RESULTS}/facets/{{comparison_id}}/segments/facets_annotated_segments.tsv", comparison_id=FACETS_COMPARISONS),
        facets_genes=expand(f"{RESULTS}/facets/{{comparison_id}}/segments/facets_annotated_genes.tsv", comparison_id=FACETS_COMPARISONS),
        qc=f"{RESULTS}/summary/qc.tsv"
    output:
        f"{RESULTS}/summary/report.html"
    log:
        f"{RESULTS}/summary/report.log"
    params:
        project=PROJECT
    conda:
        "envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {log}); "
        "python scripts/make_report.py "
        "--project {params.project} "
        "--cnvkit-segments {input.cnvkit_segments} "
        "--cnvkit-genes {input.cnvkit_genes} "
        "--facets-segments {input.facets_segments} "
        "--facets-genes {input.facets_genes} "
        "--qc {input.qc} "
        "--output {output} > {log} 2>&1"
