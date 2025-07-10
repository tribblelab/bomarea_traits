library(RevGadgets)
library(ggplot2)
library(ape)
setwd("~/Desktop/bomarea_traits/")
binary_type <- processAncStates("output/binary_type_ase_ard.tree", 
                                state_labels = c("0" = "simple",
                                                 "1" = "compound"))

plotAncStatesPie(binary_type, tip_labels = FALSE) +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22)
  )

ggsave("figures/binary_type_asr_tree.png", width = 15, height = 10, dpi = 200)
