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

g <- plotAncStatesPie(branchiness, tip_labels = TRUE)
print(g)
ggsave("figures/branchiness_2rates_asr_tree.png", width = 10, height = 10, dpi = 200)
