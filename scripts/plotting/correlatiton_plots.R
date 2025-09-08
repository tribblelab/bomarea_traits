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

