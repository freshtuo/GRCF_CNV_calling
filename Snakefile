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
        # QC is comparison-scoped because BaseSpace/DRAGEN may produce a
        # distinct BAM for the same biological control in different analyses.
        expand(
            f"{RESULTS}/qc/{{comparison_id}}/quickcheck.{{sample_id}}.txt",
            zip,
            comparison_id=QC_COMPARISON_IDS,
            sample_id=QC_SAMPLE_IDS,
        ),
        expand(
            f"{RESULTS}/qc/{{comparison_id}}/flagstat.{{sample_id}}.txt",
            zip,
            comparison_id=QC_COMPARISON_IDS,
            sample_id=QC_SAMPLE_IDS,
        ),
        expand(
            f"{RESULTS}/qc/{{comparison_id}}/idxstats.{{sample_id}}.txt",
            zip,
            comparison_id=QC_COMPARISON_IDS,
            sample_id=QC_SAMPLE_IDS,
        ),
        # CNVkit annotation runs once per comparison with run_cnvkit=yes.
        expand(CNVKIT_ANNOTATED_SEGMENTS, comparison_id=CNVKIT_COMPARISONS),
        expand(CNVKIT_ANNOTATED_GENES, comparison_id=CNVKIT_COMPARISONS),
        # FACETS runs only for comparison rows with run_facets=yes.
        expand(FACETS_SEGMENTS, comparison_id=FACETS_COMPARISONS),
        expand(FACETS_ANNOTATED_SEGMENTS, comparison_id=FACETS_COMPARISONS),
        expand(FACETS_ANNOTATED_GENES, comparison_id=FACETS_COMPARISONS),
        # Summary/report files keep project-level indexes and per-comparison pages.
        f"{RESULTS}/summary/all_comparisons.qc.tsv",
        f"{RESULTS}/summary/all_comparisons.cnv_summary.tsv",
        f"{RESULTS}/summary/all_comparisons.purity_ploidy.tsv",
        f"{RESULTS}/summary/report.html"
