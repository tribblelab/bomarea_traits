library(tidyverse)
df_list <- list(type_dat, size_dat, flowers_dat, branchiness_dat, density_dat, elevation_dat)
combined_df <- df_list %>% reduce(full_join, by = "label")

sd(combined_df$avg_size)

