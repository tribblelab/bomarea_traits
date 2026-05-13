setwd("~/Desktop/bomarea_traits/scripts/latent_liability/")

source("latent_liability_figs_functions.R")
library(RColorBrewer)

extractPrecision(
  log.file.path = "size_density/bomarea_latent_liability_size_density.log",
  out.file.path = "size_density/bomarea_latent_liability_size_density_PRECISION.log"
)

extractCorrelations(
  log.file.path = "size_density/bomarea_latent_liability_size_density.log",
  out.file.path = "size_density/bomarea_latent_liability_size_density_CORR.log",
  trait.names   = c("type", "size", "density")
)

extractVarCovar(
  log.file.path = "size_density/bomarea_latent_liability_size_density.log",
  out.file.path = "size_density/bomarea_latent_liability_size_density_VARCOVAR.log",
  trait.names = c("type", "size", "density")
)

corrs <- read.table(
  "size_density/bomarea_latent_liability_size_density_CORR.log",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

pdf("~/Desktop/bomarea_traits/figures/latent_liability_size_density.pdf", width = 7, height = 7)
plotCorrelationDistributions(
  correlations = corrs,
  labels = c("type", "size", "density"),
  lwd = 2,
  cex.labels = 1.2,
  omi = c(1,1,1,1)
)
dev.off()

varcov <- read.table(
  "bomarea_latent_liability_VARCOVAR.log",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

type_var     <- varcov$type_with_type
size_var     <- varcov$size_with_size
size_flower <- varcov$size_flower_with_size_flower

plotLatentLiabilityVarianceComparison(
  variances = list(type_var, size_var, size_flower_var),
  discrete.or.continuous = c("d", "c"),
  trait.names = c("type", "density")
)
