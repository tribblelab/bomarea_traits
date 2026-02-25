library(dplyr)
library(tidyr)
library(readxl)
library(ape)
library(car)

setwd("~/Desktop/bomarea_traits/")
source("scripts/functions.R")

traits <- read_xlsx("data_prep/bomarea_traits.xlsx", sheet = 1, na = "N/A")

## Calculates max number of branches per species
flowers <- traits %>%
  mutate(max_flowers_branch = 1 + pmax(numBranch1, numBranch2, numBranch3,
                            numBranch4, numBranch5, na.rm = TRUE),
         flowers = log(numBranchP * max_flowers_branch)) %>%
  group_by(acceptedName) %>%
  summarise(flowers = max(flowers, na.rm = TRUE)) %>%
  ungroup()

flowers$acceptedName <- gsub(" ", "_", flowers$acceptedName)

## Tree and data cleaning
# match names to tree and drop some tips
tree <- read.tree("data_prep/bom_only_MAP.tre")
tips_to_drop <- grep(
  paste0(
    "caudata|herbertiana|glaucescens|parvifolia|",
    "tribachiata|angustipetala|lehmannii|killipii|",
    "chimborazensis|trimorphophylla|hartwegii|",
    "alstroemeriodes|superba|acuminata|enanorojo|",
    "foliosa|straminea|pauciflora|distichophylla|",
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
tree_df <- left_join(tree_df, flowers,
                     by = c("speciesName" = "acceptedName"))
flowers_dat <- data.frame(label = tree_df$label,
                              flowers = tree_df$flowers)
flowers_dat <- flowers_dat[is.na(flowers_dat$label) == FALSE, ]


# manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  flowers = c(3, 6, 12)
)

# merge into existing df
flowers_dat <- flowers_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(flowers = coalesce(flowers.y, flowers.x)) %>%
  select(label, flowers)


# matrix for flowers
flowers_mat <- as.matrix(flowers_dat$flowers)
rownames(flowers_mat) <- flowers_dat$label
colnames(flowers_mat) <- "flowers"

## Write to nexus
write.nexus.data(flowers_mat, file = "data/flowers.nexus",
                 format = "standard", missing = "?")