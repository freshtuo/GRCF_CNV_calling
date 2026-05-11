#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(facets)
})

option_list <- list(
  make_option("--pileup", type = "character"),
  make_option("--comparison-id", type = "character"),
  make_option("--segments", type = "character"),
  make_option("--purity-ploidy", type = "character"),
  make_option("--plot-pdf", type = "character"),
  make_option("--plot-png", type = "character"),
  make_option("--preproc-cval", type = "integer", default = 25),
  make_option("--proc-cval", type = "integer", default = 150),
  make_option("--min-nhet", type = "integer", default = 15)
)
opt <- parse_args(OptionParser(option_list = option_list))

rcmat <- readSnpMatrix(opt$pileup)
xx <- preProcSample(rcmat, cval = opt$`preproc-cval`)
oo <- procSample(xx, cval = opt$`proc-cval`, min.nhet = opt$`min-nhet`)
fit <- emcncf(oo, min.nhet = opt$`min-nhet`)

pdf(opt$`plot-pdf`)
plotSample(x = oo, emfit = fit)
dev.off()

png(opt$`plot-png`, width = 1800, height = 1200, res = 150)
plotSample(x = oo, emfit = fit)
dev.off()

segs <- as.data.table(fit$cncf)
if (!"chrom" %in% names(segs) && "chromosome" %in% names(segs)) {
  setnames(segs, "chromosome", "chrom")
}
fwrite(segs, opt$segments, sep = "\t")

purity <- if (!is.null(fit$purity)) fit$purity else NA_real_
ploidy <- if (!is.null(fit$ploidy)) fit$ploidy else NA_real_
emflags <- if (!is.null(fit$emflags)) trimws(paste(fit$emflags, collapse = " ")) else ""
fwrite(
  data.table(
    comparison_id = opt$`comparison-id`,
    purity = purity,
    ploidy = ploidy,
    purity_emflags = emflags,
    facets_preproc_cval = opt$`preproc-cval`,
    facets_proc_cval = opt$`proc-cval`,
    facets_min_nhet = opt$`min-nhet`
  ),
  opt$`purity-ploidy`,
  sep = "\t"
)
