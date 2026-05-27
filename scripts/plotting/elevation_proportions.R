library(tidyverse)
library(ape)

setwd("~/Desktop/bomarea_traits/")

elev_nex <- read.nexus.data("data/elevation.nexus")

elev_df <- data.frame(
  label = names(elev_nex),
  elevation_state = toupper(as.character(unlist(elev_nex))),
  stringsAsFactors = FALSE
)

elev_df$species <- sapply(
  strsplit(elev_df$label, "_"),
  function(x) paste(x[1], x[2])
)

elev_df <- elev_df %>%
  filter(!grepl("^\\S+ (cf|sp)$", species, ignore.case = TRUE))

elev_lookup <- data.frame(
  elevation_state = c("0", "1", "2", "3", "4", "5",
                      "6", "7", "8", "9", "A", "B"),
  elevation_bin = c("0", "01", "012", "1", "12", "123",
                    "2", "23", "234", "3", "34", "4"),
  stringsAsFactors = FALSE
)

elev_df <- elev_df %>%
  left_join(elev_lookup, by = "elevation_state")

elev_long <- elev_df %>%
  mutate(elevation_bin = strsplit(elevation_bin, "")) %>%
  unnest(elevation_bin) %>%
  mutate(
    elevation_bin = recode(
      elevation_bin,
      "0" = "0–1000",
      "1" = "1000–2500",
      "2" = "2500–3200",
      "3" = "3200–4000",
      "4" = ">4000"
    ),
    elevation_bin = factor(
      elevation_bin,
      levels = c("0–1000",
                 "1000–2500",
                 "2500–3200",
                 "3200–4000",
                 ">4000")
    )
  )

type_nex <- read.nexus.data("data/type.nexus")

type_df <- data.frame(
  label = names(type_nex),
  type_state = as.character(unlist(type_nex)),
  stringsAsFactors = FALSE
)

type_df$species <- sapply(
  strsplit(type_df$label, "_"),
  function(x) paste(x[1], x[2])
)

type_df <- type_df %>%
  filter(!grepl("^\\S+ (cf|sp)$", species, ignore.case = TRUE)) %>%
  mutate(
    type = ifelse(type_state == "0", "Umbel", "Compound"),
    type = factor(type, levels = c("Umbel", "Compound"))
  ) %>%
  distinct(species, type)

plot_df <- elev_long %>%
  left_join(type_df, by = "species") %>%
  filter(
    !is.na(type),
    !is.na(elevation_bin)
  )


prop_df <- plot_df %>%
  distinct(species, elevation_bin, type) %>%
  count(elevation_bin, type) %>%
  group_by(elevation_bin) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


p <- ggplot(prop_df,
            aes(x = elevation_bin,
                y = proportion,
                fill = type)) +
  geom_col(
    position = "fill",
    colour = "black",
    linewidth = 0.8
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = scales::label_number(accuracy = 0.01),
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = c(
      "Compound" = "#16b4b8",
      "Umbel" = "#ef7269"
    )
  ) +
  labs(
    title = "Proportion of inflorescence types across elevation",
    x = "Elevation (m)",
    y = "Proportion of Species",
    fill = "Inflorescence Type"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain"),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(colour = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ggsave(
  "figures/elevation_type_proportions.pdf",
  plot = p,
  width = 10,
  height = 10,
  device = "pdf"
)
