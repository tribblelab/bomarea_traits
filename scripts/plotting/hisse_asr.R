setwd("~/Desktop/bomarea_traits/")
library(ggplot2)
library(RevGadgets)
library(ape)

# read in and process the ancestral states
HiSSE_file <- paste0("output/asr_hisse_binary_type.tree")
p_anc <- processAncStates(HiSSE_file,
                          state_labels = c("0A", "0B", "1A", "1B"),
                          labels_as_numbers = FALSE)

# plot the ancestral states
plot <- plotAncStatesMAP(t = p_anc, tree_layout="rectangular",
                         state_transparency = 0.5,
                         node_size = c(0.1, 5),
                         tip_labels_size = 2,
                         tip_states_size=2,
                         node_color = c("#D9081D","#F8D3D8","#072AC8","#A2D6F9")) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.1,0.85),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.title = element_text(size=6), #change legend title font size
        legend.text = element_text(size=4))

ggsave(paste0("figures/HiSSE_anc_states_activity_period.png"),plot, width=8, height=8)

# read in and process the log file
HiSSE_file <- paste0("output/hisse_tutorial.log")
pdata <- processSSE(HiSSE_file)

# plot the rates
plot <- plotMuSSE(pdata) +
  theme(legend.position = c(0.875,0.915),
        legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6))

ggsave(paste0("figures/HiSSE_div_rates_activity_period.png"),plot, width=5, height=5)

anc_states <- processAncStates(path ="output/asr_hisse_binary_type.tree",state_labels=c("0"="insect A","1"="wind A","2"="insect B", "3"="wind B"))
plotAncStatesMAP(t = anc_states, tree_layout="rectangular",
                 state_transparency = 0.5,
                 node_size = c(0.1, 5),
                 tip_labels_size = 2,
                 tip_states_size=2,
                 node_color = c("#D9081D","#F8D3D8","#072AC8","#A2D6F9"))