setwd("~/Desktop/bomarea_traits/")
library(ggplot2)
library(RevGadgets)
library(ape)

# read in and process the ancestral states
HiSSE_file <- paste0("output/asr_hisse_binary_type_timetree.tree")
p_anc <- processAncStates(HiSSE_file,
                          state_labels = c(
                            "0" = "simple",
                          # "1" = "simple_B",
                            "1" = "compound",
                            "2" = "compound_B"
                          ))

# plot the ancestral states
plot <- plotAncStatesMAP(t = p_anc, tree_layout="rectangular",
                         state_transparency = 0.75,
                         node_size = c(0.1, 5),
                         tip_labels_size = 2.3,
                         tip_states_size=3,
                         node_color = c(
                           "simple" = "#F25278",
                        #  "simple_B" = "#FFA7BA",
                           "compound" = "#7D8BE0",
                           "compound_B" = "#A2D6F9")) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.1,0.85),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6))
plot
ggsave(paste0("figures/HiSSE_asr_known_tips.png"),plot, width=8, height=8)