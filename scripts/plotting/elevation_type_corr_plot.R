library(RevGadgets)
library(ggplot2)

CHARACTER_A <- "elevation"
CHARACTER_B <- "type"

# Specify the input file
file <- paste0("output/", CHARACTER_A, "_", CHARACTER_B, "_corr_RJ.log")

# Read the trace and discard burnin
trace_qual <- readTrace(path = file, burnin = 0.25)

# Bayes factors and thresholds
BF <- c(3.2, 10, 100)
thresholds <- BF / (1 + BF)

# All variables to plot (20 total)
vars_to_plot <- c(
  "rate_A0_A1_when_B0", "rate_A0_A1_when_B1",
  "rate_A1_A0_when_B0", "rate_A1_A0_when_B1",
  "rate_A1_A2_when_B0", "rate_A1_A2_when_B1",
  "rate_A2_A1_when_B0", "rate_A2_A1_when_B1",
  "rate_A2_A3_when_B0", "rate_A2_A3_when_B1",
  "rate_A3_A2_when_B0", "rate_A3_A2_when_B1",
  "rate_B0_B1_when_A0", "rate_B0_B1_when_A1",
  "rate_B0_B1_when_A2", "rate_B0_B1_when_A3",
  "rate_B1_B0_when_A0", "rate_B1_B0_when_A1",
  "rate_B1_B0_when_A2", "rate_B1_B0_when_A3"
)

# Split variables into two groups (max 10 per group)
vars_group1 <- vars_to_plot[1:10]
vars_group2 <- vars_to_plot[11:20]

# Plot first group
plot1 <- plotTrace(trace = trace_qual, vars = vars_group1)[[1]] +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, linetype = "solid", color = "black") +
  geom_hline(yintercept = thresholds, linetype = c("longdash", "dashed", "dotted"), color = "red") +
  geom_hline(yintercept = 1 - thresholds, linetype = c("longdash", "dashed", "dotted"), color = "red") +
  theme(legend.position = c(0.40, 0.825))

# Plot second group
plot2 <- plotTrace(trace = trace_qual, vars = vars_group2)[[1]] +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, linetype = "solid", color = "black") +
  geom_hline(yintercept = thresholds, linetype = c("longdash", "dashed", "dotted"), color = "red") +
  geom_hline(yintercept = 1 - thresholds, linetype = c("longdash", "dashed", "dotted"), color = "red") +
  theme(legend.position = c(0.40, 0.825))

# Save plots
ggsave(paste0("figures/", CHARACTER_A, "_", CHARACTER_B, "_corr_RJ_part1.png"),
       plot1, width = 5, height = 5)

ggsave(paste0("figures/", CHARACTER_A, "_", CHARACTER_B, "_corr_RJ_part2.png"),
       plot2, width = 5, height = 5)
