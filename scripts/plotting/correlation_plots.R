library(ggplot2)
library(dplyr)
library(ape)

setwd("~/Desktop/bomarea_traits/")

######################
## Plot 2 varibales ##
######################

combined <- sparsity_dat %>%
  left_join(size_dat, by = "label") %>%
  left_join(flowers_dat, by = "label")


combined_subset <- combined[combined$sparsity>0,]

plot(combined_subset$sparsity,
     combined_subset$avg_size,
     xlab = "Sparsity",
     ylab = "Average Size",
     main = "Sparsity vs. Avg Size",
     pch = 16)

# cor(combined_subset$sparsity, combined_subset$avg_size, use = "complete.obs")
# cor(combined_subset$sparsity, combined_subset$avg_size, 
    # method = "spearman", use = "complete.obs")

######################
## Plot 3 varibales ##
######################

combined <- sparsity_dat %>%
  left_join(size_dat, by = "label") %>%
  left_join(flowers_dat, by = "label")%>%
  left_join(type_dat, by = "label")%>%
  left_join(sparsity_dat, by = "label")

# combined$max_branch <- as.factor(combined$max_branch)

ggplot(combined, aes(x = flowers, y = sparsity.x, color = type)) +
  geom_point(size = 3) +
  labs(x = "Number of Flowers",
       y = "Avg Size",
       color = "Inflorescence Type",
       title = "# of Flowers vs. Avg Size (inflorescence color)")



###############################################################################
## Histogram of flower number vs # of individuals (inflorescence type color) ##
###############################################################################

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




##############################################################
## Histogram Elevation vs # of species (inflorescence color ##
##############################################################

library(dplyr)
library(ape)
setwd("~/Desktop/bomarea_traits/")

merged <- left_join(elevation_dat, type_dat, by = "label")

# Label the elevation bins
merged$elevation_label <- recode_factor(
  as.character(merged$elevation_bin),
  "0" = "0–1000",
  "1" = "1000–2500",
  "2" = "2500–4500",
  "3" = "4500+"
)

# Summarize counts per elevation zone and inflorescence type
summary_df <- combined %>%
  group_by(elevation_label, type) %>%
  summarise(species_count = n())

# Bar plot
ggplot(merged, aes(x = elevation_label, fill = type)) +
  geom_bar(position = "dodge", color = "black") +
  labs(
    title = "Elevation and # of species (inflorescence type)",
    x = "Elevation (m)",
    y = "Number of Species",
    fill = "Inflorescence Type"
  ) +
  theme_classic(base_size = 14)



