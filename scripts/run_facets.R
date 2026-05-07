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
  make_option("--plot-png", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

rcmat <- readSnpMatrix(opt$pileup)
xx <- preProcSample(rcmat)
oo <- procSample(xx)
fit <- emcncf(oo)

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
fwrite(
  data.table(comparison_id = opt$`comparison-id`, purity = purity, ploidy = ploidy),
  opt$`purity-ploidy`,
  sep = "\t"
)
