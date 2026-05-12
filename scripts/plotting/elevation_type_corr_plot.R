library(RevGadgets)
library(ggplot2)
library(tidyr)
library(ggthemes)

setwd("~/Desktop/bomarea_traits/")

rates <- readTrace(c("output/corr/elevation_type_corr_RJ_run_1.log",
                     "output/corr/elevation_type_corr_RJ_run_2.log"))

rates <- combineTraces(rates)

df <- rates[[1]]

cols <- c("rate_decrease_when_0",
          "rate_decrease_when_1",
          "rate_increase_when_0",
          "rate_increase_when_1")


# change column names to more informative parameters
#    A  B
# A  X  1
# B  2  X
# where A = plain, B = compound

df_cols <- df[, cols]
colnames(df_cols) <- c("Decrease when umbellate",
                       "Decrease when compound",
                       "Increase when umbellate",
                       "Increase when compound")

df_cols %>%
  tidyr::gather(key = "grp",
                value = "val",
                factor_key = TRUE) -> df_rates


colors <- c("#D62828",
            "#F77F00",
            "#FCBF49",
            "#0085CC")

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
  stat_summary(fun = mean,
               aes(x = grp, y = val),
               geom = "point",
               size = 2,
               color = "black") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(name = "Transition type") +
  ylab("Transition rate") +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 14,
                                   angle = 30,
                                   hjust = 1),
        axis.title.y = element_text(face = "bold",
                                    size = 20),
        axis.title.x = element_text(face = "bold",
                                    size = 20))

summarizeTrace(rates,
               vars = c("rate_decrease_when_0",
                        "rate_decrease_when_1",
                        "rate_increase_when_0",
                        "rate_increase_when_1"))

print(g)
ggsave("figures/elevation_type_violin.pdf", width = 10, height = 15, dpi = 200)
