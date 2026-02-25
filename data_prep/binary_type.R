library(dplyr)
library(readxl)
library(ape)
library(tidytree)

setwd("~/Desktop/bomarea_traits/")
source("scripts/functions.R")

traits <- read_xlsx("data_prep/bomarea_traits.xlsx", sheet = 1, na = "N/A")

## Determining inflorescence type
# branching and bracteole presence
traits %>%
  group_by(acceptedName) %>%
  summarise(numBranchP = mean(numBranchP, na.rm = TRUE),
            numBracts = mean(numBracts, na.rm = TRUE),
            numBranch1 = mean(numBranch1, na.rm = TRUE),
            numBranch2 = mean(numBranch2, na.rm = TRUE),
            numBranch3 = mean(numBranch3, na.rm = TRUE),
            numBranch4 = mean(numBranch4, na.rm = TRUE),
            numBranch5 = mean(numBranch5, na.rm = TRUE),
            Bracteoles = most_frequent_discrete_value(Bracteoles)) -> trait_data

# add to trait_data
trait_data <- cbind(trait_data[, c("acceptedName",
                                   "numBranchP",
                                   "numBracts",
                                   "numBranch1")],
                    trait_data[, c("Bracteoles")])

# measures inflorescence type and cleaning species name
trait_data %>%
  mutate(umbellike = numBranch1==0) %>%
  mutate(bracteoles = Bracteoles=="Y") %>%
  mutate(acceptedName = gsub(" ", "_", acceptedName)) %>%
  select(acceptedName, umbellike, bracteoles) -> type_df

type_df$type = apply(type_df, 1, typeset)


## Tree and data cleaning
# match names to tree and drop some tips
tree <- read.tree("data/bom_only_MAP.tre")
tips_to_drop <- grep(
  paste0(
    "caudata|herbertiana|glaucescens|parvifolia|",
    "tribachiata|angustipetala|lehmannii|killipii|",
    "chimborazensis|trimorphophylla|hartwegii|",
    "alstroemeriodes|superba|acuminata|enanorojo|",
    "foliosa|straminea|pauciflora|distichophylla|",
    "Bomarea_edulis_Brazil_Campbell8900|",
    "Bomarea_edulis_Venezuela_Bunting4817" # |Bomarea_ovata_Peru_Farfan526" # don't include if doing correlated
  ),
  tree$tip.label,
  value = TRUE
)
tree_edited <- ape::drop.tip(tree, tips_to_drop)
write.tree(tree_edited, file = "data/tree_edited.tre") # change to tree_edited_corr.tre if doing correlated
tree_df <- as_tibble(tree_edited)

# combine species names
tree_df$speciesName <- unlist(lapply(tree_df$label, get_gen_sp))
tree_df <- left_join(tree_df, type_df,
                     by = c("speciesName" = "acceptedName"))
type_dat <- data.frame(label = tree_df$label,
                           type = tree_df$type)
type_dat <- type_dat[is.na(type_dat$label) == FALSE, ]

# manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  type = c(0, 2, 2)
)

# merge into existing df
type_dat <- type_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(type = ifelse(!is.na(type.y), type.y, type.x)) %>%
  select(label, type)

# turn all bracteoles into simple
# (previously, type was divided into 3 categories, this will reduce data to binary)
type_dat$type <- ifelse(type_dat$type == "1", "0", type_dat$type)
type_dat$type <- ifelse(type_dat$type == 2, 1, type_dat$type)

# make a matrix for type
type_mat <- matrix(type_dat$type, ncol = 1)
rownames(type_mat) <- type_dat$label
colnames(type_mat) <- "type"

## Write to nexus

# for everything else
write.nexus.data(type_mat, file = "data/binary_type.nexus",
                 format = "standard", missing = "?")

# for correlation analysis
# write.nexus.data(type_mat, file = "data/binary_type_corr.nexus",
#                  format = "standard", missing = "?")