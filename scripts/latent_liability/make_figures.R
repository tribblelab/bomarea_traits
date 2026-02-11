setwd("~/Desktop/bomarea_traits/scripts/latent_liability/")

source("make_latent_liability_figures(1).R")
library(RColorBrewer)

extractPrecision(
  log.file.path = "density/bomarea_latent_liability.log",
  out.file.path = "density/bomarea_latent_liability_PRECISION.log"
)

extractCorrelations(
  log.file.path = "density/bomarea_latent_liability_PRECISION.log",
  out.file.path = "density/bomarea_latent_liability_CORR.log",
  trait.names   = c("type", "density")
)

extractVarCovar(
  log.file.path = "density/bomarea_latent_liability_PRECISION.log",
  out.file.path = "density/bomarea_latent_liability_VARCOVAR.log",
  trait.names = c("type", "density")
)

corrs <- read.table(
  "density/bomarea_latent_liability_CORR.log",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

plotCorrelationDistributions(
  correlations = corrs,
  labels = c("type", "density"),
  lwd = 2,
  cex.labels = 1.2,
  omi = c(1,1,1,1)
)

varcov <- read.table(
  "bomarea_latent_liability_VARCOVAR.log",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

type_var     <- varcov$type_with_type
size_var     <- varcov$size_with_size
density_var <- varcov$density_with_density

plotLatentLiabilityVarianceComparison(
  variances = list(type_var, size_var, density_var),
  discrete.or.continuous = c("d", "c"),
  trait.names = c("type", "density")
)
