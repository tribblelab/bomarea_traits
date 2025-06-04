library(RevGadgets)
library(ggplot2)
library(ape)
setwd("~/Desktop/bomarea_traits/")
branchiness <- processAncStates("output/binary_type_ase_ard.tree", 
                                state_labels = c("0" = "simple",
                                                 "1" = "compound"))

plotAncStatesPie(infl_type, tip_labels = TRUE)
ggsave("binary_type_asr_tree.png", width = 10, height = 10, dpi = 200)