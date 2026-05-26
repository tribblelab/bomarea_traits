library(RevGadgets)

# read the output
samples <- readTrace("output/multivariate_size_sparsity_0.5.log")

# plot the posterior distribution
pdf("figures/correlations.pdf", height=4)
plotTrace(samples, vars=paste0("correlations[",1,"]"))
dev.off()
