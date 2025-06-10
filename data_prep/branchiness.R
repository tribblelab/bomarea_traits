library(dplyr)
library(tidyr)
library(readxl)
library(ape)
library(car)

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

# Calculates max number of branches per species
branchiness_max <- traits %>%
  mutate(max_branch = pmax(numBranch1, numBranch2, numBranch3,
                           numBranch4, numBranch5, na.rm = TRUE)) %>%
  group_by(acceptedName) %>%
  summarise(max_branch = max(max_branch, na.rm = TRUE)) %>%
  ungroup()

branchiness_max$acceptedName <- gsub(" ", "_", branchiness_max$acceptedName)



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
tree_df <- left_join(tree_df, branchiness_max,
                     by = c("speciesName" = "acceptedName"))
branchiness_dat <- data.frame(label = tree_df$label,
                              max_branch = tree_df$max_branch)
branchiness_dat <- branchiness_dat[is.na(branchiness_dat$label) == FALSE, ]


# Manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  max_branch = c(0, 5, 2)
)

# Merge into existing df
branchiness_dat <- branchiness_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(max_branch = coalesce(max_branch.y, max_branch.x)) %>%
  select(label, max_branch)


# Matrix for branchiness
branchiness_mat <- as.matrix(branchiness_dat$max_branch)
rownames(branchiness_mat) <- branchiness_dat$label
colnames(branchiness_mat) <- "branchiness"

# Write to nexus
write.nexus.data(branchiness_mat, file = "data/branchiness.nexus",
                 format = "standard", missing = "?")