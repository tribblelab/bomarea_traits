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


# Calculates size
size_df <- traits %>%
  mutate(totalBranchLength = lengthTotal1 +
      coalesce(lengthSeg1_1, 0) + coalesce(lengthBranch1_1, 0) +
      coalesce(lengthSeg1_2, 0) + coalesce(lengthBranch1_2, 0) +
      coalesce(lengthSeg1_3, 0) + coalesce(lengthBranch1_3, 0) +
      coalesce(lengthSeg1_4, 0) + coalesce(lengthBranch1_4, 0) +
      coalesce(lengthSeg1_5, 0),
    size = totalBranchLength / numBranchP
  ) %>%
  group_by(acceptedName) %>%
  summarize(
    avg_size = mean(size, na.rm = TRUE)
  ) %>%
  ungroup()

size_df$acceptedName <- gsub(" ", "_", size_df$acceptedName)
size_df <- size_df %>%
  mutate(size = round(avg_size, 2))

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
tree_df <- left_join(tree_df, size_df,
                     by = c("speciesName" = "acceptedName"))
size_dat <- data.frame(label = tree_df$label,
                       avg_size = tree_df$size)
size_dat <- size_dat[is.na(size_dat$label) == FALSE, ]


# Manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  size = c(5.718, 695.805, 140.205)
)

# Merge into existing df
size_dat <- size_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(size = coalesce(avg_size.y, avg_size.x)) %>%
  select(label, avg_size)


# Matrix for size
size_mat <- as.matrix(size_dat$size)
rownames(size_mat) <- size_dat$label
colnames(size_mat) <- "size"

# Write to nexus
write.nexus.data(size_mat, file = "data/siz.nexus",
                 format = "standard", missing = "?")