setwd("~/Desktop/bomarea_traits/")

library(dplyr)
library(readxl)
library(ape)
library(tidytree)

traits <- read_xlsx("data_prep/bomarea_traits.xlsx", sheet = 1, na = "N/A")

# Function to get most discrete value
most_frequent_discrete_value <- function(vec) {
  tbl <- table(vec)
  if (dim(tbl) > 0) {
    value <- names(tbl)[which.max(tbl)]
    return(value)
  } else {
    return(NA)
  }
}

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


# Calculates max branching and sparsity
traits$maxBranch <- pmax(traits$numBranch1,
                         traits$numBranch2,
                         traits$numBranch3,
                         traits$numBranch4,
                         traits$numBranch5, na.rm = TRUE)

sparsity_df <- traits %>%
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
         sparsity = ifelse(numBranchesMeasured == 0, 0,
                    (log((totalBranchLength / numBranchesMeasured) / (maxBranch + 1)))) 
  ) %>%
  select(acceptedName, sparsity) %>%
  ungroup()

sparsity_df$acceptedName <- gsub(" ", "_", sparsity_df$acceptedName)
sparsity_df <- sparsity_df %>%
  mutate(sparsity = round(sparsity, 2))

# Match names to tree and drop some tips
tree <- read.tree("data/bom_only_MAP.tre")
tips_to_drop <- grep(
  paste0(
    "caudata|herbertiana|glaucescens|parvifolia|",
    "tribachiata|angustipetala|lehmannii|",
    "chimborazensis|trimorphophylla|hartwegii|",
    "alstroemeriodes|superba|acuminata|enanorojo|",
    "killipii|foliosa|straminea|pauciflora|distichophylla|",
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
tree_df <- left_join(tree_df, sparsity_df,
                     by = c("speciesName" = "acceptedName"))
sparsity_dat <- data.frame(label = tree_df$label,
                           sparsity = tree_df$sparsity)
sparsity_dat <- sparsity_dat[is.na(sparsity_dat$label) == FALSE, ]


# Manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  sparsity = c(0.757, 2.842, 2.146)
)

# Merge into existing df
sparsity_dat <- sparsity_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(sparsity = coalesce(sparsity.y, sparsity.x)) %>%
  select(label, sparsity)


# Matrix for sparsity
sparsity_mat <- as.matrix(sparsity_dat$sparsity)
rownames(sparsity_mat) <- sparsity_dat$label
colnames(sparsity_mat) <- "sparsity"

# Write to nexus
write.nexus.data(sparsity_mat, file = "data/sparsity.nexus",
                 format = "standard", missing = "?")
