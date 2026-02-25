library(dplyr)
library(readxl)
library(ape)
library(tidytree)

setwd("~/Desktop/bomarea_traits/")
source("scripts/functions.R")

traits <- read_xlsx("data_prep/bomarea_traits.xlsx", sheet = 1, na = "N/A")

## Calculates size
size_df <- traits %>%
  mutate(numBranchesMeasured =
            rowSums(!is.na(select(., starts_with("lengthSeg"))))) %>%
  mutate(totalBranchLength = lengthTotal1 +
      coalesce(lengthSeg1_1, 0) + coalesce(lengthBranch1_1, 0) +
      coalesce(lengthSeg1_2, 0) + coalesce(lengthBranch1_2, 0) +
      coalesce(lengthSeg1_3, 0) + coalesce(lengthBranch1_3, 0) +
      coalesce(lengthSeg1_4, 0) + coalesce(lengthBranch1_4, 0) +
      coalesce(lengthSeg1_5, 0),
    size = ifelse(numBranchesMeasured == 0, 0,
           (log((totalBranchLength / numBranchesMeasured) * numBranchP)))
  ) %>%
  group_by(acceptedName) %>%
  summarize(
    avg_size = mean(size, na.rm = TRUE)
  ) %>%
  ungroup()

size_df$acceptedName <- gsub(" ", "_", size_df$acceptedName)

## Tree and data cleaning
# match names to tree and drop some tips
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
tree_df <- left_join(tree_df, size_df,
                     by = c("speciesName" = "acceptedName"))
size_dat <- data.frame(label = tree_df$label,
                       avg_size = tree_df$avg_size)
size_dat <- size_dat[is.na(size_dat$label) == FALSE, ]


# manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  avg_size = c(0.757, 2.842, 2.146)
)

# merge into existing df
size_dat <- size_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(avg_size = coalesce(avg_size.y, avg_size.x)) %>%
  select(label, avg_size)

size_dat <- size_dat %>%
        mutate(avg_size = round(avg_size, 2))

# matrix for size
size_mat <- as.matrix(size_dat$avg_size)
rownames(size_mat) <- size_dat$label
colnames(size_mat) <- "avg_size"

## Write to nexus
write.nexus.data(size_mat, file = "data/size.nexus",
                 format = "continuous", missing = "?")
