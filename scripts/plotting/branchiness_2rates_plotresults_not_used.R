library(RevGadgets)
library(ggplot2)
library(ape)
setwd("~/Desktop/bomarea_traits/")
branchiness <- processAncStates("output/branchiness_ase_ard.tree", 
                                state_labels = c("0" = "no branching",
                                                 "1" = "1",
                                                 "2" = "2",
                                                 "3" = "3",
                                                 "4" = "4",
                                                 "5" = "5",
                                                 "7" = "7"))

plotAncStatesPie(branchiness, tip_labels = FALSE) +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22)
  )

ggsave("figures/branchiness_2rates_nolabels_asr_tree.png", width = 15, height = 10, dpi = 200)
