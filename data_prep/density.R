library(dplyr)
library(readxl)
library(ape)
library(tidytree)

setwd("~/Desktop/bomarea_traits/")
source("scripts/functions.R")

traits <- read_xlsx("data_prep/bomarea_traits.xlsx", sheet = 1, na = "N/A")

# Calculates max branching and density
traits$maxBranch <- pmax(traits$numBranch1,
                         traits$numBranch2,
                         traits$numBranch3,
                         traits$numBranch4,
                         traits$numBranch5, na.rm = TRUE)

density_df <- traits %>%
  mutate(numBranchesMeasured =
           rowSums(!is.na(select(., starts_with("lengthSeg"))))) %>%
  group_by(acceptedName) %>%
  slice_max(maxBranch, with_ties = FALSE) %>%
  mutate(totalBranchLength = lengthTotal1 +
         coalesce(lengthSeg1_1, 0) +
         coalesce(lengthBranch1_1, 0) + coalesce(lengthSeg1_2, 0) +
         coalesce(lengthBranch1_2, 0) + coalesce(lengthSeg1_3, 0) +
         coalesce(lengthBranch1_3, 0) + coalesce(lengthSeg1_4, 0) +
         coalesce(lengthBranch1_4, 0) + coalesce(lengthSeg1_5, 0),
         density = ifelse(numBranchesMeasured == 0, 0,
                    (log(((totalBranchLength / numBranchesMeasured) * numBranchP) / (numBranchP * (maxBranch + 1)))) 
  )) %>%
  select(acceptedName, density) %>%
  ungroup()

density_df$acceptedName <- gsub(" ", "_", density_df$acceptedName)

## Tree and data cleaning
# Match names to tree and drop some tips
tree <- read.tree("data_prep/bom_only_MAP.tre")
tips_to_drop <- grep(
  paste0(
    "caudata|herbertiana|glaucescens|parvifolia|",
    "tribachiata|angustipetala|lehmannii|",
    "chimborazensis|trimorphophylla|hartwegii|",
    "alstroemeriodes|superba|acuminata|enanorojo|",
    "killipii|foliosa|straminea|pauciflora|distichophylla|",
    "Bomarea_edulis_Brazil_Campbell8900|",
    "Bomarea_edulis_Venezuela_Bunting4817|",
    "Bomarea_multiflora_CultivatedinCAfromCol_Greenhouse"
  ),
  tree$tip.label,
  value = TRUE
)
tree_edited <- ape::drop.tip(tree, tips_to_drop)
write.tree(tree_edited, file = "data/tree_edited.tre")
tree_df <- as_tibble(tree_edited)


# combine species names
tree_df$speciesName <- unlist(lapply(tree_df$label, get_gen_sp))
tree_df <- left_join(tree_df, density_df,
                     by = c("speciesName" = "acceptedName"))
density_dat <- data.frame(label = tree_df$label,
                           density = tree_df$density)
density_dat <- density_dat[is.na(density_dat$label) == FALSE, ]


# manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  density = c(log((((10.22 + 13.18 + 14.72) / 3) * 3) / (3 * (0 + 1))),
              log((((570.99 + 577.51 + 408.38) / 3) * 9) / (9 * (5 + 1))),
              log((((179.4 + 198.3) / 2) * 3) / (6 * (2 + 1)))
))

# merge into existing df
density_dat <- density_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(density = coalesce(density.y, density.x)) %>%
  select(label, density) %>%
  mutate(density = round(density, 2))

# matrix for density
density_mat <- as.matrix(density_dat$density)
rownames(density_mat) <- density_dat$label
colnames(density_mat) <- "density"

## Write to nexus
write.nexus.data(density_mat, file = "data/density.nexus",
                 format = "standard", missing = "?")
