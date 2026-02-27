library(RevGadgets)
library(tidyverse)
library(ggthemes)
setwd("~/Desktop/bomarea_traits/")

# Process data
rates <- readTrace(c("output/branchiness_ard_run_1.log",
                     "output/branchiness_ard_run_2.log"))

rates <- combineTraces(rates)

df <- rates[[1]]
cols <- c("er[1]", "er[2]", "er[3]",
          "er[4]", "er[5]", "er[6]",
          "er[7]")

#    A  B  C  D  E  F  G  H
# A  X  1  2  3  4  5  6  7
# B  8  X  9 10 11 12 13 14
# C 15 16  X 17 18 19 20 21
# D 22 23 24  X 25 26 27 28
# E 29 30 31 32  X 33 34 35
# F 36 37 38 39 40  X 41 42
# G 43 44 45 46 47 48  X 49
# H 50 51 52 53 54 55 56  X
# A = 0, B = 1, C = 2, D = 3, E = 4, F = 5, G = 6, H = 7

df_cols <- df[,cols]
colnames(df_cols) <- c("0 to 1",
                       "1 to 2",
                       "2 to 3",
                       "3 to 4",
                       "4 to 5",
                       "5 to 6",
                       "6 to 7")


df_cols %>%
  tidyr::gather(key = "grp",
                value = "val",
                factor_key = TRUE) -> df_rates


# set up colors
colors <- RevGadgets::colFun(7)

names(colors) <- levels(df_rates$grp)

g <- ggplot(df_rates) +
  geom_violin(data = df_rates,
              aes(x = grp,
                  y = val,
                  group = grp,
                  fill = grp),
              color = "black",
              lwd = 1,
              scale = "width",
              show.legend = FALSE) +
  stat_summary(fun = mean, aes(x = grp, y = val),
               geom = "point", size = 2, color = "black") +
  #geom_hline(yintercept = 0.0, color = "grey") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(name = "Transition rate") +
  ylab("Posterior density") +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 14,
                                   hjust = .5),
        axis.title.y = element_text(face = "bold",
                                    size = 20),
        axis.title.x = element_text(face = "bold",
                                    size = 20))

#summarizeTrace(rates, vars = c("er[1]", "er[2]", "er[3]",
#                                "er[4]", "er[5]", "er[6]",
#                                "er[7]"))

print(g)

ggsave("/figures/branchiness_ase_ard_violinPlot.png", width = 15, height = 10, dpi = 200)
