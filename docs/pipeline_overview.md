# GRCF CNV Calling Pipeline Overview

This note describes the pipeline at two levels:

1. The Snakemake workflow layer: how metadata becomes jobs.
2. The tool layer: how CNVkit, FACETS, and annotation are configured for each job.

## 1. Snakemake Workflow Layer

The workflow entry point is `Snakefile`.

```text
Snakefile
  -> config/config.yaml
  -> workflow/rules/common.smk
  -> workflow/rules/validation.smk
  -> workflow/rules/bam_qc.smk
  -> workflow/rules/cnvkit.smk
  -> workflow/rules/facets.smk
  -> workflow/rules/annotation.smk
  -> workflow/rules/report.smk
```

`common.smk` loads the project name, output directory, metadata tables, and shared helper functions.

```python
PROJECT = config["project"]
OUTDIR = config.get("outdir", "results")
RESULTS = f"{OUTDIR}/{PROJECT}"
SAMPLES_TSV = config.get("samples_tsv", "config/samples.tsv")
COMPARISONS_TSV = config.get("comparisons_tsv", "config/comparisons.tsv")
```

The two metadata tables have different roles:

- `samples.tsv` has one row per physical sample and stores BAM/BAI paths.
- `comparisons.tsv` has one row per analysis and controls which comparisons run.

`common.smk` turns those TSVs into dictionaries:

```python
SAMPLES = {row["sample_id"]: row for row in _read_tsv(SAMPLES_TSV)}
COMPARISONS = {row["comparison_id"]: row for row in _read_tsv(COMPARISONS_TSV)}
```

Then it creates the active comparison lists:

```python
CNVKIT_COMPARISONS = [
    cid for cid, row in COMPARISONS.items() if _is_yes(row.get("run_cnvkit", ""))
]

FACETS_COMPARISONS = [
    cid for cid, row in COMPARISONS.items() if _is_yes(row.get("run_facets", ""))
]
```

So `comparisons.tsv` is the main workflow driver. If a row has `run_cnvkit=yes`, CNVkit targets are created for that `comparison_id`. If a row has `run_facets=yes`, FACETS targets are created for that `comparison_id`.

## Rule All

`rule all` defines the final requested outputs. Snakemake works backward from these files to decide which rules must run.

The final targets include:

- metadata validation
- comparison-scoped BAM QC for each case/control BAM used by each comparison
- CNVkit annotations for every `run_cnvkit=yes` comparison
- FACETS segments and annotations for every `run_facets=yes` comparison
- merged project summary tables
- project-level HTML index and one HTML report per comparison

Conceptually:

```text
samples.tsv + comparisons.tsv + config.yaml
        |
        v
validate_metadata
        |
        +--> BAM QC per comparison/sample
        |
        +--> CNVkit per run_cnvkit comparison
        |       |
        |       +--> cnvkit call
        |       +--> scatter plot
        |       +--> diagram plot
        |       |
        |       v
        |   annotate CNVkit segments
        |
        +--> FACETS per run_facets comparison
                |
                v
            annotate FACETS segments

annotated CNVkit + annotated FACETS + QC
        |
        v
summarize_results
        |
        v
make_report
```

## Wildcards

Output patterns such as:

```python
CNVKIT_CNR = f"{RESULTS}/cnvkit/{{comparison_id}}/cnr/{{comparison_id}}.cnr"
```

become this after Python evaluates the f-string:

```text
results/<project>/cnvkit/{comparison_id}/cnr/{comparison_id}.cnr
```

Snakemake treats `{comparison_id}` as a wildcard. When it needs a specific file, it fills in the actual value, such as:

```text
results/test_cnv_project/cnvkit/P1_T1_vs_B1/cnr/P1_T1_vs_B1.cnr
```

Inside rule code, the value is accessed as:

```python
wildcards.comparison_id
```

That is why rules use code like:

```python
case_bam=lambda wildcards: sample_bam(comparison_case_id(wildcards.comparison_id))
```

`comparison_id` is not a normal Python variable. It is a Snakemake wildcard value known only when Snakemake creates a specific job.

## 2. Tool Layer

### Metadata Validation

`validate_metadata.py` checks the metadata before expensive jobs run.

It verifies:

- required columns exist in `samples.tsv` and `comparisons.tsv`
- `sample_id` and `comparison_id` values are unique
- `comparison_id` values are path-safe
- every sample has BAM and BAI paths
- BAM and BAI paths exist
- each comparison's `case_id` exists in `samples.tsv`
- optional `control_id` exists if provided
- `case_id` and `control_id` are not identical when a control is provided
- case/control species match the comparison species
- species resources exist in `config.yaml`
- CNVkit is globally enabled if any row has `run_cnvkit=yes`
- FACETS is only enabled for human paired non-tumor-only comparisons

The output is:

```text
results/{project}/metadata/validated.ok
```

Most downstream rules depend on this file, so failed metadata validation stops the workflow early.

### BAM QC

BAM QC runs once per comparison/sample pair. This is intentional because
BaseSpace/DRAGEN can produce distinct BAM files for the same biological control
when that control is included in multiple tumor-normal analyses.

The QC rules run:

```bash
samtools quickcheck -v sample.bam
samtools flagstat sample.bam
samtools idxstats sample.bam
```

Outputs are written under:

```text
results/{project}/qc/{comparison_id}/quickcheck.{sample_id}.txt
results/{project}/qc/{comparison_id}/flagstat.{sample_id}.txt
results/{project}/qc/{comparison_id}/idxstats.{sample_id}.txt
```

### CNVkit

CNVkit runs once per comparison with `run_cnvkit=yes`.

The rule determines the case BAM from:

```python
comparison_case_id(wildcards.comparison_id)
sample_bam(case_id)
```

If `control_id` is present, the rule adds:

```bash
--normal control.bam
```

If `control_id` is blank, the rule passes a bare `--normal`, which tells CNVkit
to build a flat reference for tumor-only analysis.

The main command is:

```bash
cnvkit.py batch case.bam [--normal control.bam | --normal] \
    --method wgs \
    --fasta reference.fa \
    --access access.bed \
    --output-dir workdir
```

The reference FASTA and access BED come from species-specific config resources:

```yaml
resources:
  human:
    fasta: ...
    access_bed: ...
```

CNVkit names output files from BAM basenames, so the rule copies the discovered `.cnr` and `.cns` files into stable comparison-based names:

```text
results/{project}/cnvkit/{comparison_id}/cnr/{comparison_id}.cnr
results/{project}/cnvkit/{comparison_id}/cns/{comparison_id}.cns
```

Then downstream CNVkit rules run:

```bash
cnvkit.py call input.cns -o output.call.cns
cnvkit.py scatter input.cnr -s input.cns -o scatter.pdf
cnvkit.py scatter input.cnr -s input.cns -o scatter.png
cnvkit.py scatter input.cnr -s input.cns -c chrN -o chromosomes/chrN.scatter.png
cnvkit.py diagram input.cns -o diagram.pdf
```

The scatter plot is kept in both PDF and PNG formats. PDF is useful for
editing, but dense WGS scatter PDFs can be slow to open because they contain
many vector points. PNG is usually easier for quick review. The workflow also
writes per-chromosome scatter PNGs for focused inspection.

### FACETS

FACETS runs once per comparison with `run_facets=yes`.

Validation expects FACETS comparisons to be:

- human
- paired
- not tumor-only
- globally enabled in `config.yaml`

The rule uses:

```python
case_bam = sample_bam(comparison_case_id(wildcards.comparison_id))
control_bam = sample_bam(comparison_control_id(wildcards.comparison_id))
common_snps = config["resources"]["human"].get("common_snps_vcf", "")
```

The first command builds a tumor/normal allele-count matrix:

```bash
snp-pileup --gzip --pseudo-snps 100 common_snps.vcf.gz \
    pileup/{comparison_id}.csv.gz \
    control.bam case.bam
```

Then the R helper runs FACETS:

```bash
Rscript scripts/run_facets.R \
    --pileup pileup/{comparison_id}.csv.gz \
    --comparison-id {comparison_id} \
    --segments facets_segments.tsv \
    --purity-ploidy facets_purity_ploidy.tsv \
    --plot-pdf facets.pdf \
    --plot-png facets.png
```

`run_facets.R` performs:

```r
rcmat <- readSnpMatrix(opt$pileup)
xx <- preProcSample(rcmat)
oo <- procSample(xx)
fit <- emcncf(oo)
```

It writes:

```text
results/{project}/facets/{comparison_id}/segments/facets_segments.tsv
results/{project}/facets/{comparison_id}/purity_ploidy/facets_purity_ploidy.tsv
results/{project}/facets/{comparison_id}/plots/facets.pdf
results/{project}/facets/{comparison_id}/plots/facets.png
```

### Annotation

Annotation is handled by `scripts/annotate_segments.py`.

CNVkit annotation uses:

```bash
python scripts/annotate_segments.py \
    --source-tool cnvkit \
    --segments {comparison_id}.call.cns \
    --genes-bed species.genes.bed \
    --thresholds gain,loss,amp,deep_del,focal_mb,broad_mb \
    --out-segments annotated_segments.tsv \
    --out-genes annotated_genes.tsv
```

FACETS annotation uses the same script with:

```bash
--source-tool facets
--segments facets_segments.tsv
```

The script normalizes segment columns from CNVkit or FACETS into a common schema:

```text
chromosome
segment_start
segment_end
segment_length
log2_ratio
total_cn
minor_cn
cnv_call
loh_status
event_size
```

It also harmonizes chromosome naming before the gene overlap step. For example,
FACETS may emit chromosomes as `1`, while the gene BED may use `chr1`; the
annotation script makes those styles compatible before running `bedtools
intersect`.

For CNVkit, calls are based on log2 thresholds from `config.yaml`:

```yaml
cnv_thresholds:
  gain: 0.3
  loss: -0.3
  amp: 1.0
  deep_del: -1.0
```

For FACETS, calls are based on total copy number:

```text
total_cn >= 5  -> amplification
total_cn >= 3  -> gain
total_cn <= 0  -> homozygous_deletion
total_cn <= 1  -> loss
```

LOH is assigned when:

```text
minor_cn == 0
```

Event size is assigned using:

```yaml
event_size:
  focal_mb: 5
  broad_mb: 10
```

The script overlaps each segment with `genes_bed` and writes:

```text
cnvkit/{comparison_id}/annotation/annotated_segments.tsv
cnvkit/{comparison_id}/annotation/annotated_genes.tsv
facets/{comparison_id}/annotation/annotated_segments.tsv
facets/{comparison_id}/annotation/annotated_genes.tsv
```

## Summary And Report

`summarize_results` merges:

- CNVkit annotated segment/gene counts
- FACETS annotated segment/gene counts
- FACETS purity/ploidy when available
- comparison-scoped QC outputs

It writes:

```text
results/{project}/summary/all_comparisons.qc.tsv
results/{project}/summary/all_comparisons.cnv_summary.tsv
results/{project}/summary/all_comparisons.purity_ploidy.tsv
```

`make_report` creates:

```text
results/{project}/summary/report.html
results/{project}/summary/reports/{comparison_id}.report.html
```

`summary/report.html` is a project-level index. Each comparison also gets its
own detailed report under `summary/reports/`.

## Dry-Run Check

A dry run with the current example metadata builds a 58-job DAG:

```text
job                    count
-------------------  -------
all                        1
annotate_cnvkit            6
annotate_facets            3
cnvkit_batch               6
cnvkit_call                6
cnvkit_diagram             6
cnvkit_scatter             6
facets                     3
make_report                1
samtools_flagstat          6
samtools_idxstats          6
samtools_quickcheck        6
summarize_results          1
validate_metadata          1
total                     58
```

This means the workflow structure is coherent enough for Snakemake to plan all jobs.

However, the current checked-in config and metadata contain placeholder paths such as:

```text
/path/hg38.fa
/path/blood.bam
```

A real run will fail metadata validation until those are replaced with real paths.

## DAG Plot Commands

The usual Snakemake commands are:

```bash
conda run -n snakemake snakemake --dag --cores 1 > docs/workflow_dag.dot
dot -Tpdf docs/workflow_dag.dot -o docs/workflow_dag.pdf
dot -Tpng docs/workflow_dag.dot -o docs/workflow_dag.png
```

For a simpler rule-level graph:

```bash
conda run -n snakemake snakemake --rulegraph --cores 1 > docs/workflow_rulegraph.dot
dot -Tpdf docs/workflow_rulegraph.dot -o docs/workflow_rulegraph.pdf
dot -Tpng docs/workflow_rulegraph.dot -o docs/workflow_rulegraph.png
```

In this environment, `snakemake -n -p --cores 1` works, but Snakemake 8.2.3 raises an internal `execution_settings` error for `--dag` and `--rulegraph`. Because of that, `docs/workflow_rulegraph.dot` was written manually from the parsed rule dependencies and rendered with Graphviz.

## Things To Double-Check

Before a real run:

- Replace all placeholder BAM/BAI/resource paths.
- Confirm `resources.human.common_snps_vcf` exists if any comparison has `run_facets=yes`.
- Confirm CNVkit `--method wgs` is appropriate for the BAMs. For targeted/panel/WES data, this may need a different CNVkit setup.
- Review FACETS defaults in `scripts/run_facets.R`; currently `preProcSample`, `procSample`, and `emcncf` use package defaults.
- Confirm gene BED chromosome naming is compatible with the references. The annotation script can harmonize common `chr1` versus `1` differences.
- Interpret tumor-only CNVkit results cautiously because germline CNVs cannot be subtracted without a matched control.
