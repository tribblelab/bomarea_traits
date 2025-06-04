library(RevGadgets)
library(ggplot2)
library(ape)
setwd("~/Desktop/bomarea_traits/")
branchiness <- processAncStates("output/binary_type_ase_ard.tree", 
                                state_labels = c("0" = "simple",
                                                 "1" = "compound"))

plotAncStatesPie(infl_type, tip_labels = TRUE)
ggsave("binary_type_asr_tree.png", width = 10, height = 10, dpi = 200)

#rates <- readTrace(c("output/binary_type_ard_run_1.log",
#                     "output/binary_type_ard_run_2.log"))

#rates <- combineTraces(rates)

#plotTrace(rates, match = "rate")

library(RevGadgets)
setwd("~/Desktop/bomarea_infl_testruns/")
infl_type <- processAncStates("output/infl_type_ase_ard.tree", state_labels = c("0" = "Plain umbel",
                                                                                "1" = "Bracteole umbel",
                                                                                "2" = "Compound inflorescence"))
plotAncStatesPie(infl_type, tip_labels = FALSE)