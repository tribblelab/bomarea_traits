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

elevation_min_max <- elevation_min_max %>%
  rename(label = acceptedName)

type_dat <- type_dat %>%
  mutate(label = sub("^(Bomarea_[A-Za-z]+).*", "\\1", label))

merged <- left_join(elevation_min_max, type_dat, by = "label")

# Label the elevation bins
merged$elevation_label <- recode_factor(
  as.character(merged$elevation_bin),
  "0" = "0–1000",
  "1" = "1000–2500",
  "2" = "2500–3200",
  "3" = "3200-4500",
  "4" = "4500-5000"
)

# Summarize counts per elevation bin and inflorescence type
summary_df <- merged %>%
  group_by(, type) %>%
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


library(dplyr)
library(tidyr)
library(ggplot2)

# --- Step 1. Define elevation bins ---
bins <- data.frame(
  elevation_label = c("0–1000", "1000–2500", "2500–4000", "4000-5000"),
  lower = c(0, 1000, 2500, 4000),
  upper = c(1000, 2500, 4000, 5000)
)

# --- Step 2. Expand each species into overlapping bins ---
expanded <- tidyr::crossing(collapsed, bins) %>%
  filter(max_elevation > lower & min_elevation < upper) %>%
  select(label, elevation_label, type)

# --- Step 3. Remove species with missing type values ---
expanded <- expanded %>%
  filter(!is.na(type))

# --- Step 4. Summarize counts per elevation bin and type ---
summary_df <- expanded %>%
  group_by(elevation_label, type) %>%
  summarise(species_count = n_distinct(label), .groups = "drop")

# --- Step 5. Convert counts to ratios (sum = 1.00 for each elevation bin) ---
summary_df <- summary_df %>%
  group_by(elevation_label) %>%
  mutate(ratio = species_count / sum(species_count)) %>%
  ungroup()

# --- Step 6. Plot as stacked histogram with proportions (0–1 scale) ---
ggplot(summary_df, aes(x = elevation_label, y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(
    title = "Proportion of inflorescence types across elevation",
    x = "Elevation (m)",
    y = "Proportion of Species",
    fill = "Inflorescence Type"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )
