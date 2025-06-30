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
traits$elevation_avg <- (traits$min_elevation + traits$max_elevation) / 2
traits$elevation_bin <- cut(traits$elevation_avg,
                        breaks = c(-Inf, 1000, 2500, 4500, Inf),
                        labels = c("0", "1", "2", "3"),
                        right = FALSE)
elevation_df <- traits[, c("acceptedName", "elevation_avg", "elevation_bin")]
elevation_df$acceptedName <- gsub(" ", "_", elevation_df$acceptedName)

elevation_df <- elevation_df %>%
  group_by(acceptedName) %>%
  summarise(
    elevation_avg = first(elevation_avg),
    elevation_bin = first(elevation_bin)
  )

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
    "Bomarea_edulis_Venezuela_Bunting4817"
  ),
  tree$tip.label,
  value = TRUE
)
tree_edited <- ape::drop.tip(tree, tips_to_drop)
write.tree(tree_edited, file = "data/tree_edited.tre")
tree_df <- as_tibble(tree_edited)


# Combine species names
tree_df$speciesName <- unlist(lapply(tree_df$label, get_gen_sp))
tree_df <- left_join(tree_df, elevation_df,
                     by = c("speciesName" = "acceptedName"))
elevation_dat <- data.frame(label = tree_df$label,
                            elevation_bin = tree_df$elevation_bin)
elevation_dat <- elevation_dat[is.na(elevation_dat$label) == FALSE, ]


# Manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  elevation_bin = c(2, 2, 2)
)

# Merge into existing df
elevation_dat <- elevation_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(elevation_bin = coalesce(as.character(elevation_bin.y), as.character(elevation_bin.x))) %>%
  select(label, elevation_bin)


# Matrix for elevation
elevation_mat <- as.matrix(elevation_dat$elevation_bin)
rownames(elevation_mat) <- elevation_dat$label
colnames(elevation_mat) <- "elevation"

# Write to nexus
write.nexus.data(elevation_mat, file = "data/elevation.nexus",
                 format = "standard", missing = "?")
