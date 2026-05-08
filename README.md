# GRCF_CNV_calling

Snakemake workflow for copy-number analysis from existing DRAGEN/BaseSpace BAMs.
The workflow does not realign reads. CNVkit is the baseline caller for all
enabled comparisons; FACETS is optional for human paired comparisons.

## Layout

```text
Snakefile
config/
  config.yaml
  samples.tsv
  comparisons.tsv
workflow/rules/
  validation.smk
  bam_qc.smk
  cnvkit.smk
  facets.smk
  annotation.smk
  report.smk
envs/
scripts/
resources/
results/
```

## Metadata

`config/samples.tsv` has one row per physical sample and records BAM/BAI paths.
`sample_type` is descriptive only and does not control workflow logic.

`config/comparisons.tsv` has one row per CNV analysis and drives all execution:

- `case_id`: sample being evaluated
- `control_id`: matched baseline/reference sample, blank for tumor-only
- `run_cnvkit`: run CNVkit for this comparison
- `run_facets`: run FACETS for this comparison
- `comparison_type`: descriptive reporting label

Validation checks that case/control samples exist, species are consistent, BAM
and BAI paths exist, and FACETS is only enabled for human paired non-tumor-only
comparisons.

`comparison_id` values are used in output paths, so they should be unique and
path-safe: letters, numbers, dots, dashes, and underscores only. Validation also
rejects duplicate sample/comparison IDs and comparisons where `case_id` and
`control_id` are the same.

## Configure

Edit `config/config.yaml` and replace the placeholder paths:

- `resources.<species>.fasta`
- `resources.<species>.access_bed`
- `resources.<species>.genes_bed`
- `resources.human.common_snps_vcf` if any comparison has `run_facets=yes`

For WGS CNVkit runs, `cnvkit.target_avg_size` controls the bin size passed to
`--target-avg-size`. Set it to blank/null to use CNVkit's automatic bin sizing.

Then replace example BAM/BAI paths in `config/samples.tsv`.

## Run

```bash
snakemake -n --cores 1
snakemake --use-conda --cores 8
```

Outputs are written under:

```text
results/{project}/
  qc/{comparison_id}/
  cnvkit/{comparison_id}/
    annotation/
  facets/{comparison_id}/
    annotation/
  summary/
    reports/
```

The summary outputs are:

- `summary/all_comparisons.qc.tsv`: comparison-scoped BAM QC
- `summary/all_comparisons.cnv_summary.tsv`: CNV event counts by comparison and caller
- `summary/all_comparisons.purity_ploidy.tsv`: FACETS purity/ploidy when available
- `summary/report.html`: project-level index
- `summary/reports/{comparison_id}.report.html`: one detailed report per comparison

Tumor-only CNV calls should be interpreted cautiously because germline CNVs
cannot be removed without a matched control.
