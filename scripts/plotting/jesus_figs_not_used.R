library(tidyverse)
df_list <- list(type_dat, size_dat, flowers_dat, branchiness_dat, density_dat, elevation_dat)
combined_df <- df_list %>% reduce(full_join, by = "label")


library(tidyverse)
library(ape)
library(ggtree)

setwd("~/Desktop/bomarea_traits/")

tree <- read.tree("data/tree_edited.tre")

traits <- combined_df %>%
  filter(!is.na(type)) %>%
  mutate(
    type = factor(
      type,
      levels = c(0, 1),
      labels = c("umbel", "compound")
    ),
    unlogged_flowers = exp(flowers)
  )

# continuous traits
cont_traits <- c("avg_size", "unlogged_flowers", "max_branch", "density")

# keep only taxa shared between tree and data
traits <- traits %>%
  filter(label %in% tree$tip.label)

tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, traits$label)
)

# output folder
dir.create(
  "figures/jesus",
  recursive = TRUE,
  showWarnings = FALSE
)

# histos

for (tr in cont_traits) {
  
  p <- ggplot(traits, aes(x = .data[[tr]])) +
    geom_histogram(
      bins = 15,
      fill = "grey70",
      color = "black"
    ) +
    theme_bw() +
    labs(
      x = tr,
      y = "count",
      title = paste("histogram of", tr)
    )
  
  ggsave(
    paste0("figures/jesus/hist_", tr, ".pdf"),
    p,
    width = 5,
    height = 4
  )
}

# type barplot

p_type <- ggplot(traits, aes(x = type)) +
  geom_bar(
    fill = "grey70",
    color = "black"
  ) +
  theme_bw() +
  labs(
    x = "inflorescence type",
    y = "count"
  )

ggsave(
  "figures/jesus/type_barplot.pdf",
  p_type,
  width = 5,
  height = 4
)

# continuous traits on tree

for (tr in cont_traits) {
  
  plot_df <- traits %>%
    select(label, value = all_of(tr))
  
  p <- ggtree(tree) %<+% plot_df +
    geom_tippoint(
      aes(color = value),
      size = 2
    ) +
    scale_color_gradient(
      low = "lightgrey",
      high = "black",
      na.value = "grey80"
    ) +
    theme_tree2() +
    labs(
      color = tr,
      title = paste(tr, "mapped on tree")
    )
  
  ggsave(
    paste0("figures/jesus/tree_", tr, ".pdf"),
    p,
    width = 8,
    height = 10
  )
}

# type on tree

p_tree_type <- ggtree(tree) %<+% traits +
  geom_tippoint(
    aes(color = type),
    size = 2
  ) +
  theme_tree2() +
  labs(color = "inflorescence type")

ggsave(
  "figures/jesus/tree_type.pdf",
  p_tree_type,
  width = 8,
  height = 10
)

# pairwise scatterplots

trait_pairs <- combn(cont_traits, 2, simplify = FALSE)

for (pair in trait_pairs) {
  
  xvar <- pair[1]
  yvar <- pair[2]
  
  p <- ggplot(
    traits,
    aes(
      x = .data[[xvar]],
      y = .data[[yvar]],
      color = type
    )
  ) +
    geom_point(size = 2, alpha = 0.8) +
    theme_bw() +
    labs(
      x = xvar,
      y = yvar,
      color = "type",
      title = paste(xvar, "vs", yvar)
    )
  
  ggsave(
    paste0(
      "figures/jesus/",
      xvar,
      "_vs_",
      yvar,
      ".pdf"
    ),
    p,
    width = 5,
    height = 5
  )
}

# traits by inflorescence type

traits_long <- traits %>%
  select(
    type,
    avg_size,
    unlogged_flowers,
    max_branch,
    density
  ) %>%
  pivot_longer(
    cols = c(avg_size, unlogged_flowers, max_branch, density),
    names_to = "trait",
    values_to = "value"
  )

p_by_type <- ggplot(
  traits_long,
  aes(x = type, y = value, color = type)
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.7,
    size = 1.5
  ) +
  facet_wrap(
    ~trait,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    x = "inflorescence type",
    y = "trait value"
  )

ggsave(
  "figures/jesus/traits_by_type.pdf",
  p_by_type,
  width = 9,
  height = 6
)

