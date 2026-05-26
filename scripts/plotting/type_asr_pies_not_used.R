library(RevGadgets)
library(ggplot2)
library(ape)
setwd("~/Desktop/bomarea_traits/")
type <- processAncStates("output/type_ase_ard.tree", 
                                state_labels = c("0" = "simple",
                                                 "1" = "compound"))

plotAncStatesPie(type, tip_labels = FALSE) +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22)
  )

ggsave("figures/type_asr_tree.png", width = 15, height = 10, dpi = 200)
