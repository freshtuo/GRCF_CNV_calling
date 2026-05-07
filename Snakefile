configfile: "config/config.yaml"

# The Snakefile is the workflow entry point. It loads config.yaml, includes the
# modular rule files, and defines the final files requested by `rule all`.
include: "workflow/rules/common.smk"
include: "workflow/rules/validation.smk"
include: "workflow/rules/bam_qc.smk"
include: "workflow/rules/cnvkit.smk"
include: "workflow/rules/facets.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/maintenance.smk"


# `rule all` is the default target when running `snakemake`.
# Snakemake works backward from these expected outputs to decide which rules run.
rule all:
    input:
        # Metadata validation is a prerequisite for every analysis step.
        f"{RESULTS}/metadata/validated.ok",
        # QC runs once per physical sample listed in samples.tsv.
        expand(f"{RESULTS}/qc/samples/{{sample_id}}/quickcheck.txt", sample_id=SAMPLES.keys()),
        expand(f"{RESULTS}/qc/samples/{{sample_id}}/flagstat.txt", sample_id=SAMPLES.keys()),
        expand(f"{RESULTS}/qc/samples/{{sample_id}}/idxstats.txt", sample_id=SAMPLES.keys()),
        # CNVkit annotation runs once per comparison with run_cnvkit=yes.
        expand(f"{RESULTS}/annotation/{{comparison_id}}/annotated_segments.tsv", comparison_id=CNVKIT_COMPARISONS),
        expand(f"{RESULTS}/annotation/{{comparison_id}}/annotated_genes.tsv", comparison_id=CNVKIT_COMPARISONS),
        # FACETS runs only for comparison rows with run_facets=yes.
        expand(f"{RESULTS}/facets/{{comparison_id}}/segments/facets_segments.tsv", comparison_id=FACETS_COMPARISONS),
        # Summary/report files keep CNV caller outputs separate and summarize QC.
        f"{RESULTS}/summary/qc.tsv",
        f"{RESULTS}/summary/report.html"
