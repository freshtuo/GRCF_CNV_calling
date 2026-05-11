# CNV Pipeline Result Files

This note summarizes the main output files, what each file contains, and which
files are most useful for interpretation.

## Recommended Files To Review

For most users, start with these files under:

```text
results/<project>/
```

### Project Summary

```text
summary/report.html
```

Project-level HTML index with one row per comparison. It links to per-comparison
HTML reports and merged summary tables.

```text
summary/reports/<comparison_id>.report.html
```

Per-comparison HTML report with QC, purity/ploidy when available, CNV summary,
ranked high-overlap gene highlights, and links to detailed caller outputs.
Ranked gene highlights are shown separately for CNVkit and FACETS. By default,
each caller section shows up to 20 genes per category for gains,
losses/deletions, and LOH, requiring `gene_overlap_fraction >= 0.8`.
These values are controlled by `report.top_genes_per_category` and
`report.min_gene_overlap_fraction` in `config/config.yaml`.
Chromosome Y genes are excluded from these ranked highlights to reduce
sex-chromosome artifacts in samples without chrY. Low-priority non-coding or
predicted gene prefixes (`LOC`, `MIR`, `LINC`, `SNOR`, `SNORD`, `SNORA`, `RNU`,
`RNA5S`) are also hidden from ranked highlights. These are display filters only;
all rows remain in the detailed annotation TSVs.
If FACETS purity was not estimated, the FACETS report section includes a caution
note because allele-specific CN and LOH highlights are less reliable without a
confident purity/ploidy fit.

Note: HTML links are relative to the result folder layout. Links should keep
working if the whole `results/<project>/` directory is copied while preserving
relative paths. Links may break if only selected HTML files are copied.

### Main Interpretation Tables

```text
cnvkit/<comparison_id>/annotation/annotated_segments.tsv
facets/<comparison_id>/annotation/annotated_segments.tsv
```

Primary segment-level CNV/LOH annotation tables. These are the best files for
copy-number interpretation because CNV and LOH are called at segment level.

```text
cnvkit/<comparison_id>/annotation/annotated_genes.tsv
facets/<comparison_id>/annotation/annotated_genes.tsv
```

Gene-indexed view of segment calls. A gene row means the gene overlaps a segment
and inherits that segment's CNV/LOH call. Neutral genes are included, so use
`cnv_call` and `loh_status` to distinguish altered genes from copy-neutral
overlaps. This file is useful for searching genes of interest, but it is not an
independent gene-resolution call.

Summary fields named `n_affected_genes` count only genes with non-neutral
`cnv_call` or `loh_status == LOH`; neutral rows in `annotated_genes.tsv` are not
counted as affected.

## CNVkit Outputs

CNVkit runs in two main steps:

```text
cnvkit.py batch
cnvkit.py call
```

### Bin-Level Copy Ratio

```text
cnvkit/<comparison_id>/cnr/<comparison_id>.cnr
```

CNVkit `.cnr` file. This is the bin-level copy-ratio table.

Common columns:

```text
chromosome
start
end
gene
depth
log2
weight
```

Use this file for detailed signal/QC review and plotting. It has many rows, one
per genomic bin.

### Segmented Copy Ratio

```text
cnvkit/<comparison_id>/cns/<comparison_id>.cns
```

CNVkit `.cns` file from `cnvkit.py batch`. This is the segmented copy-ratio
table before integer copy-number calling.

Common columns:

```text
chromosome
start
end
gene
log2
depth
probes
weight
ci_lo
ci_hi
```

### Called Segments

```text
cnvkit/<comparison_id>/calls/<comparison_id>.call.cns
```

Output from `cnvkit.py call`. This keeps the segmented copy-ratio information
and adds integer copy number, usually in the `cn` column.

This is the CNVkit source file used to create:

```text
cnvkit/<comparison_id>/annotation/annotated_segments.tsv
cnvkit/<comparison_id>/annotation/annotated_genes.tsv
```

For interpretation, the annotated tables are usually easier to use than raw
`.call.cns`.

### CNVkit Native Batch Outputs

```text
cnvkit/<comparison_id>/batch_outputs/
```

Preserved CNVkit-native output files, including files such as:

```text
reference.cnn
*.targetcoverage.cnn
*.antitargetcoverage.cnn
*.cnr
*.cns
*.call.cns
*.bintest.cns
hg38.access.target.bed
hg38.access.antitarget.bed
```

These are useful for debugging or advanced review. Most users should focus on
the normalized pipeline outputs in `cnr/`, `cns/`, `calls/`, `annotation/`, and
`plots/`.

## FACETS Outputs

FACETS runs only for eligible paired human tumor-normal comparisons.

### FACETS Segment Table

```text
facets/<comparison_id>/segments/facets_segments.tsv
```

Main FACETS segment output from `fit$cncf`. FACETS directly reports
allele-specific segment calls, so there is no separate CNVkit-style `.cnr` and
`.cns` pair.

Common columns:

```text
chrom
start
end
num.mark
nhet
cnlr.median
mafR
tcn.em
lcn.em
cf.em
```

Important fields:

```text
cnlr.median
```

Copy-number log-ratio-like signal.

```text
tcn.em
```

Estimated total copy number.

```text
lcn.em
```

Estimated lesser/minor copy number. In this pipeline, `lcn.em == 0` is annotated
as LOH.

```text
cf.em
```

Estimated cellular fraction for the segment. This is FACETS' model-based
estimate of the fraction of tumor cells carrying the fitted copy-number event.
Values near `1` suggest a more clonal event; lower values suggest a possible
subclonal, weaker, or less certain event. Treat `cf.em` as supporting evidence,
not an exact biological percentage.

This file is the FACETS source file used to create:

```text
facets/<comparison_id>/annotation/annotated_segments.tsv
facets/<comparison_id>/annotation/annotated_genes.tsv
```

### FACETS Purity And Ploidy

```text
facets/<comparison_id>/purity_ploidy/facets_purity_ploidy.tsv
summary/all_comparisons.purity_ploidy.tsv
```

Contains FACETS global purity and ploidy estimates when available. The table
also records FACETS purity-estimation warning text and the FACETS model
parameters used for that run:

```text
purity_emflags
facets_preproc_cval
facets_proc_cval
facets_min_nhet
```

The default values in `config/config.yaml` match FACETS 0.6.2 defaults:

```text
preProcSample(..., cval = 25)
procSample(..., cval = 150, min.nhet = 15)
emcncf(..., min.nhet = 15)
```

These parameters are recorded so results remain interpretable if FACETS
defaults or project settings change later. Do not change them for routine runs
unless you understand the effect and keep parameter-specific sensitivity
results separate from primary pipeline outputs.

Blank purity with ploidy `2` can occur when FACETS produces segment calls but
does not return a confident global purity estimate. In that case, FACETS
segments can still be reviewed, but allele-specific calls and LOH should be
interpreted cautiously and cross-checked with plots and CNVkit.

### FACETS Diagnostic Plot

```text
facets/<comparison_id>/plots/facets.png
facets/<comparison_id>/plots/facets.pdf
```

The FACETS plot is generated by FACETS' standard `plotSample()` function. The
x-axis is genome position by chromosome.

Top panel, `log-ratio`:

```text
cnlr / total copy-number signal
```

Values near `0` are roughly copy-neutral. Positive shifts suggest gain or
amplification; negative shifts suggest loss or deletion.

Second panel, `log-odds-ratio`:

```text
allelic imbalance signal at informative SNPs
```

Values near `0` indicate balanced alleles. Shifts away from `0` suggest
allelic imbalance and can support LOH, especially when the fitted minor copy
number is `0`. A region with flat log-ratio but shifted log-odds-ratio can be
consistent with copy-neutral LOH.

Third panel, `copy number (em)`:

```text
tcn.em and lcn.em model-fitted copy numbers
```

Black points represent fitted total copy number (`tcn.em`). Red/pink points
represent fitted lesser/minor copy number (`lcn.em`). Common interpretations:

```text
tcn.em = 2, lcn.em = 1   copy-neutral diploid-like
tcn.em = 1, lcn.em = 0   single-copy deletion with LOH
tcn.em = 2, lcn.em = 0   copy-neutral LOH
tcn.em > 2               gain/amplification
lcn.em = 0               LOH in this pipeline's annotation
```

Bottom strip, `cf-em`:

```text
cf.em model-fitted cellular fraction
```

The colored vertical marks summarize segment-level `cf.em`. The color scale is
mainly a visual cue: pale/white indicates lower cellular fraction, blue
indicates increasing cellular fraction, and light orange/bisque indicates the
high end of the scale, often near-clonal values. Orange does not mean purity
failed; it means high fitted `cf.em` for those segments. Use the `cf.em` column
in `facets_segments.tsv` for exact values rather than estimating from color.

FACETS purity estimation usually benefits from clear variation in the first two
panels: copy-number shifts in log-ratio and/or allelic imbalance in
log-odds-ratio. Very flat panels may indicate a copy-number quiet tumor, low
tumor purity, or too little signal for a confident purity/ploidy fit. FACETS can
also fail to estimate purity when the signal is noisy or complex rather than
flat. If global purity is blank, interpret `tcn.em`, `lcn.em`, `cf.em`, and LOH
calls cautiously and cross-check with CNVkit and the diagnostic plots.

## Annotated Segment Columns

The shared annotated segment tables contain normalized columns from CNVkit or
FACETS:

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
n_genes
genes
max_gene_overlap_fraction
```

Key fields:

```text
log2_ratio
```

Segment copy-ratio signal. For CNVkit this comes from `log2`; for FACETS it
comes from `cnlr.median`.

```text
total_cn
```

Integer total copy number. For CNVkit this comes from `.call.cns` column `cn`.
For FACETS this comes from `tcn.em`.

```text
minor_cn
```

Minor/lesser allele copy number. For FACETS this comes from `lcn.em`. CNVkit
usually does not provide allele-specific minor copy number in this workflow.

```text
cnv_call
```

Pipeline classification of the segment.

For CNVkit, calls are based on `log2_ratio` thresholds:

```text
amplification: log2_ratio > amp threshold
gain:          log2_ratio > gain threshold
deep_deletion: log2_ratio < deep_del threshold
loss:          log2_ratio < loss threshold
neutral:       otherwise
unknown:       missing log2_ratio
```

For FACETS, calls are based on total copy number:

```text
amplification: total_cn >= 5
gain:          total_cn >= 3
homozygous_deletion: total_cn <= 0
loss:          total_cn <= 1
neutral:       otherwise
unknown:       missing total_cn
```

```text
loh_status
```

Loss of heterozygosity status:

```text
LOH:      minor_cn == 0
retained: otherwise
```

For FACETS, this is meaningful because FACETS estimates allele-specific copy
number. For CNVkit, `minor_cn` is usually unavailable, so `retained` mostly means
LOH was not estimated by CNVkit in this workflow.

```text
event_size
```

Segment size class based on configured focal and broad thresholds.

```text
n_genes, genes, max_gene_overlap_fraction
```

Collapsed gene overlap summary for the segment.

## Annotated Gene Columns

The gene-level tables contain one row per retained gene-segment overlap after
annotation filtering. A gene can appear more than once if it overlaps multiple
segments.

`annotated_genes.tsv` can have different row counts across comparisons. "All
genes" means all genes that overlap retained segments for that comparison, not
every gene in the reference genome exactly once. Row counts differ because each
comparison can have different segment boundaries, omitted/problematic regions,
caller-specific segmentation, and duplicate rows for genes that overlap multiple
segments. For a unique-gene count, count distinct `gene_name` values rather than
total rows.

Common columns:

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
gene_name
gene_start
gene_end
gene_length_bp
gene_overlap_bp
gene_overlap_fraction
gene_overlap_type
```

Interpretation:

```text
gene_overlap_bp
```

Number of bases where the gene overlaps the segment.

```text
gene_overlap_fraction
```

Fraction of the gene body overlapped by the segment:

```text
gene_overlap_bp / gene_length_bp
```

```text
gene_overlap_type
```

`complete` if the segment covers the full gene body, otherwise `partial`.

Important caveat: gene-level calls are inherited from segment-level calls. A
gene-level row should be read as:

```text
This gene overlaps a segment with this CNV/LOH status.
```

not as:

```text
This gene was independently called at gene-level resolution.
```

For high-confidence gene interpretation, prioritize genes with high
`gene_overlap_fraction`, complete overlap, strong CNV signal, and agreement
between callers when available.

## Annotation Sensitivity

Current overlap settings are:

```yaml
annotation:
  min_gene_overlap_bp: 1
  min_gene_overlap_fraction: 0.0
```

This is a maximum-sensitivity setting. Any gene with at least 1 bp overlap with
a segment can be annotated.

This is useful for discovery and gene lookup, but boundary genes should be
interpreted carefully. For publication-style high-confidence filtering, consider
post-filtering by:

```text
gene_overlap_fraction >= 0.5
gene_overlap_fraction >= 0.8
gene_overlap_type == complete
```

## Which Caller Should I Use?

For paired tumor-normal comparisons:

```text
CNVkit: robust total copy-ratio evidence.
FACETS: allele-specific copy number and LOH, especially when purity/ploidy are available.
```

When FACETS purity is unavailable, FACETS segments can still be reviewed, but
allele-specific calls should be interpreted cautiously and cross-checked with
CNVkit and plots.

For tumor-only comparisons:

```text
CNVkit only.
```

Tumor-only CNVkit calls should be interpreted cautiously because germline CNVs
cannot be subtracted without a matched normal.

## Plots

CNVkit plots:

```text
cnvkit/<comparison_id>/plots/<comparison_id>.scatter.png
cnvkit/<comparison_id>/plots/<comparison_id>.scatter.pdf
cnvkit/<comparison_id>/plots/<comparison_id>.diagram.pdf
cnvkit/<comparison_id>/plots/chromosomes/
```

FACETS plots:

```text
facets/<comparison_id>/plots/facets.png
facets/<comparison_id>/plots/facets.pdf
```

Use plots to confirm major events, noisy samples, ambiguous FACETS fits, and
candidate genes of interest.

## Practical Review Order

Recommended review flow:

1. Open `summary/report.html`.
2. Check per-comparison reports in `summary/reports/`.
3. Review `summary/all_comparisons.cnv_summary.tsv`.
4. Review `summary/all_comparisons.purity_ploidy.tsv` for FACETS confidence.
5. Use `annotation/annotated_segments.tsv` as the primary evidence table.
6. Use `annotation/annotated_genes.tsv` to search and prioritize genes.
7. Confirm important calls with CNVkit/FACETS plots.

## Rerunning Part Of The Workflow

Use `--rerun-triggers mtime` for targeted reruns when you only want Snakemake to
consider file timestamps. This avoids rerunning upstream CNV callers just
because rule code or parameters changed during report/annotation development.

Preview any targeted rerun first:

```bash
snakemake make_report --use-conda --cores 8 --rerun-triggers mtime -n -p
```

Regenerate only the HTML reports after report-layout or display-only changes:

```bash
snakemake make_report --use-conda --cores 8 --rerun-triggers mtime -R make_report
```

Regenerate merged summary tables and HTML reports after changes to
`scripts/summarize_cnv.py` or report inputs:

```bash
snakemake summarize_results make_report --use-conda --cores 8 --rerun-triggers mtime -R summarize_results make_report
```

Regenerate CNV annotations, summaries, and reports after changes to
`scripts/annotate_segments.py` or gene-overlap annotation settings:

```bash
snakemake annotate_cnvkit annotate_facets summarize_results make_report --use-conda --cores 8 --rerun-triggers mtime -R annotate_cnvkit annotate_facets summarize_results make_report
```

Do not force CNVkit or FACETS caller rules unless caller inputs or caller
parameters changed and you intentionally want to recalculate copy-number calls.

## Maintenance And Cleanup

Maintenance targets are optional and are not part of the default workflow.

Standard cleanup:

```bash
snakemake standard_cleanup --use-conda --cores 8 --rerun-triggers mtime
```

This removes disposable intermediates that are not used by downstream
annotation or reports:

```text
facets/<comparison_id>/pileup/
cnvkit/<comparison_id>/batch_tmp.*
```

FACETS pileups are regenerated if the FACETS rule is rerun. CNVkit `batch_tmp.*`
directories are temporary work directories and are only expected to remain after
failed or interrupted jobs.

Standard cleanup also gzips large CNVkit-native files inside:

```text
cnvkit/<comparison_id>/batch_outputs/
```

Compression uses `pigz` when it is available for parallel gzip-compatible
compression, controlled by `maintenance.gzip_threads` in `config/config.yaml`
(default: 4 threads per comparison). If `pigz` is not available, the workflow
falls back to standard single-threaded `gzip`.

The canonical user-facing files in `cnr/`, `cns/`, `calls/`, `annotation/`, and
`plots/` are left unchanged.

Advanced lower-level targets are also available:

```bash
snakemake cleanup_intermediates --use-conda --cores 8 --rerun-triggers mtime
snakemake gzip_cnvkit_batch_outputs --use-conda --cores 8 --rerun-triggers mtime
```

Deep cleanup:

```bash
snakemake deep_cleanup --use-conda --cores 8 --rerun-triggers mtime
```

This runs the disposable-file cleanup and also removes:

```text
cnvkit/<comparison_id>/batch_outputs/
```

Use deep cleanup only after final review, when CNVkit-native intermediate/debug
files such as `reference.cnn`, `*.targetcoverage.cnn`, native `*.cnr`, native
`*.cns`, and `*.bintest.cns` are no longer needed.
