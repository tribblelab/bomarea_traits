library(RevGadgets)
library(tidyverse)

setwd("~/Desktop/bomarea_traits/")

# files 
HiSSE_tree_file <- "output/hisse/constrained/hisse_rj_asr.tree"
HiSSE_rates_file <- "output/hisse/constrained/hisse_rj.log"
HiSSE_maps_file <- "output/hisse/constrained/hisse_rj_stoch_maps.log"

# read in data
p_tree <- readTrees(HiSSE_tree_file)
p_anc <- processAncStates(HiSSE_tree_file, state_labels = c("0" = "0A", 
                                                       "1" = "1A", 
                                                       "2" = "0B", 
                                                       "3" = "1B"))
p_rates <- readTrace(HiSSE_rates_file)  
p_maps <- processStochMaps(tree = p_tree, 
                           paths = HiSSE_maps_file,
                           states = as.character(0:3) )                

# plot the ancestral states
pies <- plotAncStatesPie(
  p_anc,
  tip_labels = TRUE,
  state_transparency = 1.0,
  node_pie_size = 0.8,
  tip_pie_size = 0.70,
  pie_colors = c(
    "0A" = "#660000",
    "1A" = "#003049",
    "0B" = "#e86064",
    "1B" = "#669BBC"
  )
) +
  theme(
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 15)
  )
pies
ggsave("~/Desktop/bomarea_traits/figures/hisse_asr.pdf",
           height = 10, width = 13)

# plot net div rates
p_rates[[1]]$net_div_0 <- p_rates[[1]]$`speciation_observed[1]` -  
                            p_rates[[1]]$`extinction_observed[1]`
p_rates[[1]]$net_div_1 <- p_rates[[1]]$`speciation_observed[2]` -  
                            p_rates[[1]]$`extinction_observed[2]`

p_rates[[1]] %>%
    select(net_div_0, net_div_1) %>%
    gather() %>%
    ggplot(aes(key, value, fill = key)) +
    geom_violin() +
    stat_summary(fun.y = "median", geom = "point", shape = 19, size = 2) +
    ylim(-0.25,1.5) + 
    scale_fill_manual(values = c(
        "net_div_0" = "#e86064",
        "net_div_1" = "#669BBC")) +
    scale_x_discrete(labels = c("Simple", "Compound")) +
    labs(x = "Net Diversification Rate", y = "Density") +
    theme_bw() +
    theme(legend.position = "none" )

ggsave("~/Desktop/bomarea_traits/figures/hisse_violin.pdf",
           height = 6, width = 6)

#plot sctochastic maps
map_plot <- plotStochMaps(tree = p_tree, 
              maps = p_maps,
              colors = c("0" = colFun(4)[1],
                         "1" = colFun(4)[3],
                         "2" = colFun(4)[2],
                         "3" = colFun(4)[4]),
              color_by = "PROBS")

ggsave("~/Desktop/bomarea_traits/figures/hisse_stoch_maps_PROBS.pdf",
           height = 10, width = 8)
