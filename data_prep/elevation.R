library(dplyr)
library(tidyr)
library(readxl)
library(ape)
library(car)

setwd("~/Desktop/bomarea_traits/")

traits <- read_xlsx("data_prep/bomarea_traits.xlsx", sheet = 1, na = "N/A")

# Function to create species name
get_gen_sp <- function(x) {
  if (is.na(x)) {
    return(NA)
  } else if (grepl("_cf_", x)) {
    namesplit <- unlist(strsplit(x, split = "_"))
    newname <- paste0(namesplit[1], "_", namesplit[3])
    return(newname)
  } else {
    namesplit <- unlist(strsplit(x, split = "_"))
    newname <- paste0(namesplit[1], "_", namesplit[2])
    return(newname)
  }
}

# Elevation bins
bins <- data.frame(
  bin_label = c("0", "1", "2", "3", "4"),
  lower = c(-Inf, 1000, 2500, 3200, 4000),
  upper = c(1000, 2500, 3200, 4000, Inf)
)

elevation_df <- crossing(traits, bins) %>%
  filter(max_elevation > lower & min_elevation < upper) %>%
  mutate(elevation_bin = bin_label) %>%
  group_by(acceptedName) %>%
  summarise(
    min_elevation = min(min_elevation, na.rm = TRUE),
    max_elevation = max(max_elevation, na.rm = TRUE),
    bins = paste0(sort(unique(elevation_bin)), collapse = ""),
    .groups = "drop"
  ) %>%
  mutate(acceptedName = gsub(" ", "_", acceptedName))

state_map <- c(
  "0"   = 0,
  "01"  = 1,
  "012" = 2,
  "1"   = 3,
  "12"  = 4,
  "123" = 5,
  "2"   = 6,
  "23"  = 7,
  "234" = 8,
  "3"   = 9,
  "34"  = "A",
  "4"   = "B"
)

elevation_df <- elevation_df %>%
  mutate(
    elevation_bin = bins,
    elevation_state = state_map[elevation_bin]
  ) %>%
  select(acceptedName, min_elevation, max_elevation, elevation_bin, elevation_state)



# traits$elevation_avg <- (traits$min_elevation + traits$max_elevation) / 2

# traits$elevation_bin <- cut(traits$elevation_avg,
                        # breaks = c(-Inf, 1000, 2500, 3200, 4500, Inf),
                        # labels = c("0", "1", "2", "3", "4"),
                        # right = FALSE)
# elevation_df <- traits[, c("acceptedName", "elevation_avg", "elevation_bin")]
# elevation_df$acceptedName <- gsub(" ", "_", elevation_df$acceptedName)

# elevation_df <- elevation_df %>%
#   group_by(acceptedName) %>%
#   summarise(
#     elevation_avg = first(elevation_avg),
#     elevation_bin = first(elevation_bin)
#   )

# Match names to tree and drop some tips
tree <- read.tree("data/bom_only_MAP.tre")
tips_to_drop <- grep(
  paste0(
    "caudata|herbertiana|glaucescens|parvifolia|",
    "tribachiata|angustipetala|lehmannii|killipii|",
    "chimborazensis|trimorphophylla|hartwegii|",
    "alstroemeriodes|superba|acuminata|enanorojo|",
    "foliosa|straminea|pauciflora|distichophylla|",
    "Bomarea_edulis_Brazil_Campbell8900|",
    "Bomarea_edulis_Venezuela_Bunting4817|Bomarea_ovata_Peru_Farfan526"
  ),
  tree$tip.label,
  value = TRUE
)
tree_edited <- ape::drop.tip(tree, tips_to_drop)
write.tree(tree_edited, file = "data/tree_edited.tre")
tree_df <- tibble(label = tree_edited$tip.label)


# Combine species names
tree_df$speciesName <- unlist(lapply(tree_df$label, get_gen_sp))
tree_df <- left_join(tree_df, elevation_df,
                     by = c("speciesName" = "acceptedName"))
elevation_dat <- data.frame(label = tree_df$label,
                            elevation_state = tree_df$elevation_state)
elevation_dat <- elevation_dat[is.na(elevation_dat$label) == FALSE, ]


# Manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  elevation_state = c(4, 4, 4)
)

# Merge into existing df
elevation_dat <- elevation_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(elevation_state = coalesce(as.character(elevation_state.y), as.character(elevation_state.x))) %>%
  select(label, elevation_state)


# Matrix for elevation
elevation_mat <- as.matrix(elevation_dat$elevation_state)
rownames(elevation_mat) <- elevation_dat$label
colnames(elevation_mat) <- "elevation"

# Write to nexus
write.nexus.data(elevation_mat, file = "data/elevation.nexus",
                 format = "standard", missing = "?")
