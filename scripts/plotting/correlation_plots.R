library(ggplot2)

# Just size and sparsity
combined <- full_join(sparsity_df,size_df)


combined_subset <- combined[combined$sparsity>0,]

plot(combined_subset$sparsity,
     combined_subset$avg_size,
     xlab = "Sparsity",
     ylab = "Average Size",
     main = "Sparsity vs. Avg Size",
     pch = 16)

cor(combined_subset$sparsity, combined_subset$avg_size, use = "complete.obs")
cor(combined_subset$sparsity, combined_subset$avg_size, 
    method = "spearman", use = "complete.obs")

# Color code (elevation,branchiness,type)

library(dplyr)

combined <- sparsity_df %>%
  full_join(size_df, by = "acceptedName") %>%
  full_join(elevation_df, by = "acceptedName")%>%
  full_join(type_df, by = "acceptedName")%>%
  full_join(branchiness_max, by = "acceptedName")

combined$max_branch <- as.factor(combined$max_branch)

ggplot(combined, aes(x = sparsity, y = avg_size, color = max_branch)) +
  geom_point(size = 3) +
  labs(x = "Sparsity",
       y = "Avg Size",
       color = "Branchiness",
       title = "Sparsity vs. Avg Size (branchiness color)")

# Overlapping histogram

library(ggplot2)
library(dplyr)


# Keep only types 0, 1, 2 and combine 1 & 2 into 1
combined <- full_join(flowers,type_df)
combined <- combined %>%
  filter(type %in% c(0,1,2)) %>%   # keep relevant types
  mutate(type_combined = case_when(
    type %in% c(0,1) ~ 0,          # type 0 & 1 → 0
    type == 2 ~ 1                   # type 2 → 1
  ),
  type_combined = factor(type_combined, levels = c(0,1))  # type 2 drawn first
  )


# Plot with gray overlay
ggplot(combined, aes(x = flowers, fill = type_combined)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.75, position = "identity") +
  scale_fill_manual(values = c("1" = "skyblue", "0" = "#E41A1C")) +
  labs(
    x = "Flower number",
    y = "Individuals",
    fill = "Inflorescence type",
    title = "Flower number and inflorescence type"
  ) +
  theme_minimal(base_size = 14)


