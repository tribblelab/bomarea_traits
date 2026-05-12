library(RevGadgets)
library(ggplot2)

CHARACTER_A <- "elevation"
CHARACTER_B <- "type"

# specify the input file
file <- paste0("output/corr/", CHARACTER_A, "_", CHARACTER_B, "_corr_RJ.log")

# read the trace and discard burnin
trace_qual <- readTrace(path = file, burnin = 0.25)

BF <- c(3.2, 10, 100)
p <- BF / (1 + BF)
# produce the plot object, showing the posterior distributions of the rates.
p <- plotTrace(
  trace = trace_qual,
  vars = c("prob_gain_A_indep", "prob_gain_B_indep", "prob_loss_A_indep", "prob_loss_B_indep")
)[[1]] +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, linetype = "solid", color = "black") +
  geom_hline(yintercept = p, linetype = c("longdash", "dashed", "dotted"), color = "red") +
  geom_hline(yintercept = 1 - p, linetype = c("longdash", "dashed", "dotted"), color = "red") +
  # modify legend location using ggplot2
  theme(legend.position = c(0.40, 0.825))

ggsave(paste0("Primates_", CHARACTER_A, "_", CHARACTER_B, "_corr_RJ.pdf"), p, width = 5, height = 5)



# confidence intervals
trace_df <- trace_qual[[1]]

apply(trace_df[, c(
  "rate_decrease_when_0",
  "rate_decrease_when_1",
  "rate_increase_when_0",
  "rate_increase_when_1"
)], 2, quantile, probs = c(0.025, 0.5, 0.975))



