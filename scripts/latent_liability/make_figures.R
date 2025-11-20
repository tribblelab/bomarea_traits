source("make_latent_liability_figures(1).R")
library(RColorBrewer)

extractPrecision(
  log.file.path = "bomarea_latent_liability.log",
  out.file.path = "bomarea_latent_liability_PRECISION.log"
)

extractCorrelations(
  log.file.path = "bomarea_latent_liability_PRECISION.log",
  out.file.path = "bomarea_latent_liability_CORR.log",
  trait.names   = c("type", "size", "sparsity")
)

extractVarCovar(
  log.file.path = "bomarea_latent_liability_PRECISION.log",
  out.file.path = "bomarea_latent_liability_VARCOVAR.log",
  trait.names = c("type", "size", "sparsity")
)

corrs <- read.table(
  "bomarea_latent_liability_CORR.log",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

plotCorrelationDistributions(
  correlations = corrs,
  labels = c("type", "size", "sparsity"),
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
sparsity_var <- varcov$sparsity_with_sparsity

plotLatentLiabilityVarianceComparison(
  variances = list(type_var, size_var, sparsity_var),
  discrete.or.continuous = c("d", "c", "c"),
  trait.names = c("type", "size", "sparsity")
)
