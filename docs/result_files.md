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
affected-gene preview, and links to detailed caller outputs.

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

Gene-indexed view of segment calls. A gene row means the gene overlaps a CNV/LOH
segment and inherits that segment's call. This is useful for searching genes of
interest, but it is not an independent gene-resolution call.

Current behavior: the gene-level file reports affected genes only, excluding
neutral and unknown segments. A planned update is to include all overlapping
genes and use `cnv_call` to distinguish neutral from altered genes.

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

Estimated cellular fraction for the segment.

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

Contains FACETS global purity and ploidy estimates when available.

Blank purity with ploidy `2` can occur when FACETS produces segment calls but
does not return a confident global purity estimate. In that case, FACETS
segments can still be reviewed, but allele-specific calls and LOH should be
interpreted cautiously and cross-checked with plots and CNVkit.

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

The gene-level tables contain one row per reported gene-segment overlap after
annotation filtering.

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

