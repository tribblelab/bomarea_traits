setwd("~/Desktop/bomarea_traits/")
devtools::install_github("revbayes/RevGadgets@stochastic_map",force=TRUE)

library(ggplot2)
library(RevGadgets)
library(ape)

setwd("~/Desktop/bomarea_traits/")

# read in and process the ancestral states
HiSSE_file <- paste0("output/hisse/hisse_rj_asr.tree")
p_anc <- processAncStates(HiSSE_file,
                          state_labels = c(
                            "0" = "simple_A",
                            "1" = "compound_A",
                            "2" = "simple_B",
                            "3" = "compound_B"))

state_colors <- c(
  "simple_A"   = "#F25278",
  "simple_B"   = "#FFA7BA",
  "compound_A" = "#7D8BE0",
  "compound_B" = "#A2D6F9"
)

# plot pies
pies <- plotAncStatesPie(
  p_anc,
  tip_labels = TRUE,
  pie_colors = state_colors,
  state_transparency = 1.0,
  node_pie_size = 1,
  tip_pie_size = 0.85
) +
  theme(
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 15)
  )

pies
ggsave("figures/hisse/HiSSE_asr_pies_test.png", pies, width = 8, height = 8)

# plot the ancestral states
plot <- plotAncStatesMAP(t = p_anc, tree_layout="rectangular",
                         state_transparency = 0.85,
                         node_size = c(0.1, 5),
                         tip_labels_size = 2.3,
                         tip_states_size=3,
                         node_color = state_colors) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.1,0.85),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6))
plot
ggsave(paste0("figures/HiSSE_asr.png"),plot, width=8, height=8)

# read in and process the log file
HiSSE_file <- paste0("output/hisse_test_run_1.log")
pdata <- processSSE(HiSSE_file)
pdata$state_label <- paste0(
  ifelse(pdata$observed_state == 0, "simple_", "compound_"),
  pdata$hidden_state
)
pdata$state <- pdata$state_label

# plot the rates
plot <- plotMuSSE(pdata) +
  aes(fill = state) +
  scale_fill_manual(values = state_colors) +
  coord_cartesian(xlim = c(-1, 2.5)) +
  theme(legend.position = c(0.875,0.915),
        legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6))

plot
ggsave(paste0("figures/hisse_div_rates_hidden_states.png"),plot, width=6, height=5)
