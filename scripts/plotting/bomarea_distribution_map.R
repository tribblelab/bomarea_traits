# ── Bomarea distribution map ──────────────────────────────────────────────────
# Reads a GBIF occurrence CSV and a NEXUS file, filters to valid species
# (excluding cf. and sp. taxa), and plots a distribution map.
#
# Required packages: ggplot2, sf, rnaturalearth, rnaturalearthdata, dplyr, ape

library(ape)        # read.nexus.data
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

setwd("~/Desktop/bomarea_traits/")

nex <- read.nexus.data("data/type.nexus")
taxa_names <- names(nex)

# Extract genus + epithet from underscore-delimited taxon labels
# e.g. "Bomarea_brevis_Bolivia_Teran1759" -> "Bomarea brevis"
species_raw <- sapply(strsplit(taxa_names, "_"), function(x) paste(x[1], x[2]))

# Exclude cf. and sp. taxa
species_clean <- unique(species_raw[!grepl("^\\S+ (cf|sp)$", species_raw,
                                           ignore.case = TRUE)])

cat("Species to map:", length(species_clean), "\n")

# ── 2. Load GBIF occurrence data ───────────────────────────────────────────────

occ <- read.delim("data_prep/occurence_data.csv", sep = "\t", quote = "",
                  stringsAsFactors = FALSE)

# Filter to target species with valid coordinates in the Americas
occ_filtered <- occ %>%
  filter(
    species %in% species_clean,
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    decimalLongitude >= -125, decimalLongitude <= -30,
    decimalLatitude  >= -55,  decimalLatitude  <= 25
  )

cat("Occurrences retained:", nrow(occ_filtered), "\n")
cat("Species with records:", n_distinct(occ_filtered$species), "\n")

# ── 3. Get base map ────────────────────────────────────────────────────────────

world <- ne_countries(scale = "medium", returnclass = "sf")

# ── 4. Plot ───────────────────────────────────────────────────────────────────

p <- ggplot() +
  # Ocean background
  theme(panel.background = element_rect(fill = "white")) +
  # Land
  geom_sf(data = world,
          fill = "#d0d0d0", colour = "#888888", linewidth = 0.3) +
  # Occurrence points
  geom_point(data = occ_filtered,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "black", alpha = 0.6, size = 0.8, shape = 16) +
  # Crop to Americas
  coord_sf(xlim = c(-110, -32), ylim = c(-50, 27), expand = FALSE) +
  # Reference lines
  geom_hline(yintercept = 0,     colour = "#888888", linewidth = 0.4,
             linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = -23.5, colour = "#888888", linewidth = 0.3,
             linetype = "dotted", alpha = 0.5) +
  annotate("text", x = -83, y = 0.7,
           label = "Equator", size = 2.2, colour = "#666666",
           hjust = 0) +
  annotate("text", x = -83, y = -22.8,
           label = "Tropic of Capricorn", size = 1.9, colour = "#666666",
           hjust = 0) +
  # Labels
  labs(
    x = "Longitude", y = "Latitude"
  ) +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(colour = "#cccccc", linewidth = 0.2),
    plot.title    = element_text(colour = "black", size = 16, face = "bold",
                                 family = "serif", hjust = 0.5,
                                 margin = margin(t = 10, b = 2)),
    plot.subtitle = element_text(colour = "#555555", size = 7.5, hjust = 0.5,
                                 margin = margin(b = 6)),
    plot.caption  = element_text(colour = "#888888", size = 5.5, hjust = 1,
                                 margin = margin(t = 4, b = 4)),
    axis.text  = element_text(colour = "#444444", size = 7),
    axis.title = element_text(colour = "#444444", size = 8),
    axis.ticks = element_line(colour = "#888888"),
    plot.margin = margin(10, 10, 10, 10)
  )

# ── 5. Save ───────────────────────────────────────────────────────────────────

ggsave("figures/bomarea_distribution.pdf", plot = p, width = 8, height = 10,
       device = "pdf")
ggsave("bomarea_distribution.png", plot = p, width = 8, height = 10,
       dpi = 250)

cat("Saved: bomarea_distribution.pdf and bomarea_distribution.png\n")
