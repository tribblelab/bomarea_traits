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


# Function to quantify branching and bracteoles
typeset <- function(df) {
  if (any(is.na(df))) {
    return(NA)
  } else {
    if (df[2]==TRUE & df[3]==FALSE) { ## umbellike, no bracteoles
      return(0)
    }
    else if (df[2]==TRUE & df[3]==TRUE) { ## umbellike w/ bracteoles
      return(1)
    }
else if (df[2]==FALSE & df[3]==TRUE) { ## non umbel (branching) w/ bracteoles
      return(2)
    } else {
      return(NA)
    }
  }
}


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

trait_data <- cbind(trait_data[, c("acceptedName",
                                   "numBranchP",
                                   "numBracts",
                                   "numBranch1")],
                    trait_data[, c("Bracteoles")])

# Measures inflorescence type
trait_data %>%
  mutate(umbellike = numBranch1==0) %>%
  mutate(bracteoles = Bracteoles=="Y") %>%
  mutate(acceptedName = gsub(" ", "_", acceptedName)) %>%
  select(acceptedName, umbellike, bracteoles) -> type_df

type_df$type = apply(type_df, 1, typeset)


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
    "Bomarea_edulis_Venezuela_Bunting4817|Bomarea_ovata_Peru_Farfan526| " #<- if doing different than corr then do include ovata
  ),
  tree$tip.label,
  value = TRUE
)
tree_edited <- ape::drop.tip(tree, tips_to_drop)
write.tree(tree_edited, file = "data/tree_edited.tre")
tree_df <- as_tibble(tree_edited)

# Combine species names
tree_df$speciesName <- unlist(lapply(tree_df$label, get_gen_sp))
tree_df <- left_join(tree_df, type_df,
                     by = c("speciesName" = "acceptedName"))
type_dat <- data.frame(label = tree_df$label,
                           type = tree_df$type)
type_dat <- type_dat[is.na(type_dat$label) == FALSE, ]

# Manually added these to nexus file
manual_add <- data.frame(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  type = c(0, 2, 2)
)

#merge into existing df
type_dat <- type_dat %>%
  left_join(manual_add, by = "label") %>%
  mutate(type = ifelse(!is.na(type.y), type.y, type.x)) %>%
  select(label, type)

# Turn all bracteoles into simple (remove for 3 types)
type_dat$type <- ifelse(type_dat$type == "1", "0", type_dat$type)

# Turn simple into 0 and compound into 1 (remove for 3 types)
type_dat$type <- ifelse(type_dat$type == 2, 1, type_dat$type)


# Make a matrix for type
type_mat <- matrix(type_dat$type, ncol = 1)
rownames(type_mat) <- type_dat$label
colnames(type_mat) <- "type"

# Write to nexus
write.nexus.data(type_mat, file = "data/binary_type_corr.nexus",
                 format = "standard", missing = "?")

