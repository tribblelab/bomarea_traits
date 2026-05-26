setwd("~/Desktop/bomarea_traits/")

library(dplyr)
library(tidyr)
library(readxl)
library(ape)
library(tidytree)

## functions
traits <- read_xlsx("data_prep/bomarea_traits.xlsx", sheet = 1, na = "N/A")

most_frequent_discrete_value <- function(vec) {
  tbl <- table(vec)
  if (length(tbl) == 0) return(NA)
  names(tbl)[which.max(tbl)]
}

get_gen_sp <- function(x) {
  if (is.na(x)) return(NA)
  if (grepl("_cf_", x)) {
    parts <- unlist(strsplit(x, "_"))
    paste0(parts[1], "_", parts[3])
  } else {
    parts <- unlist(strsplit(x, "_"))
    paste0(parts[1], "_", parts[2])
  }
}

typeset <- function(df_row) {
  if (any(is.na(df_row))) return(NA)
  umbellike  <- df_row[2]
  bracteoles <- df_row[3]
  if (umbellike == TRUE  && bracteoles == FALSE) return(0) # umbel, no bracteoles
  if (umbellike == TRUE  && bracteoles == TRUE)  return(1) # umbel, bracteoles
  if (umbellike == FALSE && bracteoles == TRUE)  return(2) # branched, bracteoles
  NA
}

## edit the tree
tree <- read.tree("data/bom_only_MAP.tre")

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

tip_tbl <- tibble(
  label = tree_edited$tip.label,
  speciesName = vapply(tree_edited$tip.label, get_gen_sp, character(1))
)

# standardize acceptedName to underscores
traits <- traits %>%
  mutate(acceptedName = gsub(" ", "_", acceptedName))

## size
size_df <- traits %>%
  mutate(numBranchesMeasured = rowSums(!is.na(select(., starts_with("lengthSeg"))))) %>%
  mutate(
    totalBranchLength = lengthTotal1 +
      coalesce(lengthSeg1_1, 0) + coalesce(lengthBranch1_1, 0) +
      coalesce(lengthSeg1_2, 0) + coalesce(lengthBranch1_2, 0) +
      coalesce(lengthSeg1_3, 0) + coalesce(lengthBranch1_3, 0) +
      coalesce(lengthSeg1_4, 0) + coalesce(lengthBranch1_4, 0) +
      coalesce(lengthSeg1_5, 0),
    size = ifelse(numBranchesMeasured == 0, 0,
                  log((totalBranchLength / numBranchesMeasured) * numBranchP))
  ) %>%
  group_by(acceptedName) %>%
  summarise(avg_size = mean(size, na.rm = TRUE), .groups = "drop") %>%
  mutate(avg_size = round(avg_size, 2))

## flowers
flowers_df <- traits %>%
  mutate(
    max_flowers_branch = 1 + pmax(numBranch1, numBranch2, numBranch3, numBranch4, numBranch5, na.rm = TRUE),
    flowers = log(numBranchP * max_flowers_branch)
  ) %>%
  group_by(acceptedName) %>%
  summarise(flowers = max(flowers, na.rm = TRUE), .groups = "drop")

## branchiness
branchiness_df <- traits %>%
  mutate(max_branch = pmax(numBranch1, numBranch2, numBranch3, numBranch4, numBranch5, na.rm = TRUE)) %>%
  group_by(acceptedName) %>%
  summarise(branchiness = max(max_branch, na.rm = TRUE), .groups = "drop")

## density
density_df <- traits %>%
  mutate(maxBranch = pmax(numBranch1, numBranch2, numBranch3, numBranch4, numBranch5, na.rm = TRUE)) %>%
  mutate(numBranchesMeasured = rowSums(!is.na(select(., starts_with("lengthSeg"))))) %>%
  group_by(acceptedName) %>%
  slice_max(maxBranch, with_ties = FALSE) %>%
  mutate(
    totalBranchLength = lengthTotal1 +
      coalesce(lengthSeg1_1, 0) + coalesce(lengthBranch1_1, 0) +
      coalesce(lengthSeg1_2, 0) + coalesce(lengthBranch1_2, 0) +
      coalesce(lengthSeg1_3, 0) + coalesce(lengthBranch1_3, 0) +
      coalesce(lengthSeg1_4, 0) + coalesce(lengthBranch1_4, 0) +
      coalesce(lengthSeg1_5, 0),
    density = ifelse(numBranchesMeasured == 0, 0,
                     log(((totalBranchLength / numBranchesMeasured) * numBranchP) /
                           (numBranchP * (maxBranch + 1))))
  ) %>%
  ungroup() %>%
  select(acceptedName, density) %>%
  mutate(density = round(density, 2))

## elevation
bins <- tibble(
  bin_label = c("0", "1", "2", "3", "4"),
  lower = c(-Inf, 1000, 2500, 3200, 4000),
  upper = c(1000, 2500, 3200, 4000, Inf)
)

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

elevation_df <- crossing(traits %>% select(acceptedName, min_elevation, max_elevation), bins) %>%
  filter(max_elevation > lower & min_elevation < upper) %>%
  mutate(elevation_bin = bin_label) %>%
  group_by(acceptedName) %>%
  summarise(
    min_elevation = min(min_elevation, na.rm = TRUE),
    max_elevation = max(max_elevation, na.rm = TRUE),
    elevation_bin = paste0(sort(unique(elevation_bin)), collapse = ""),
    .groups = "drop"
  ) %>%
  select(acceptedName, elevation_bin)

## binary type
type_df <- traits %>%
  group_by(acceptedName) %>%
  summarise(
    numBranchP = mean(numBranchP, na.rm = TRUE),
    numBracts  = mean(numBracts,  na.rm = TRUE),
    numBranch1 = mean(numBranch1, na.rm = TRUE),
    Bracteoles = most_frequent_discrete_value(Bracteoles),
    .groups = "drop"
  ) %>%
  mutate(
    umbellike  = numBranch1 == 0,
    bracteoles = Bracteoles == "Y"
  ) %>%
  select(acceptedName, umbellike, bracteoles) %>%
  mutate(type_raw = apply(., 1, typeset)) %>%
  mutate(
    type_raw = ifelse(type_raw == 1, 0, type_raw),  # bracteolate umbel -> simple
    type     = ifelse(type_raw == 2, 1, type_raw)   # branched -> compound (binary)
  ) %>%
  select(acceptedName, type)

## manual add ins
manual_add <- tibble(
  label = c("Bomarea_sp__oso_Peru_Graham12613",
            "Bomarea_sp__ponillalsoya_Peru_Graham12616",
            "Bomarea_sp__catanatasoya_Peru_Graham12611"),
  avg_size = c(0.757, 2.842, 2.146),
  flowers  = c(3, 6, 12),
  elevation_bin = c(2, 2, 2),
  density  = c(
    log((((10.22 + 13.18 + 14.72) / 3) * 3) / (3 * (0 + 1))),
    log((((570.99 + 577.51 + 408.38) / 3) * 9) / (9 * (5 + 1))),
    log((((179.4 + 198.3) / 2) * 3) / (6 * (2 + 1)))
  ),
  branchiness = c(0, 5, 2),
  type = c(0, 1, 1)
)

## master table
trait_table <- tip_tbl %>%
  left_join(type_df,       by = c("speciesName" = "acceptedName")) %>%
  left_join(branchiness_df,by = c("speciesName" = "acceptedName")) %>%
  left_join(size_df,       by = c("speciesName" = "acceptedName")) %>%
  left_join(flowers_df,    by = c("speciesName" = "acceptedName")) %>%
  left_join(density_df,    by = c("speciesName" = "acceptedName")) %>%
  left_join(elevation_df,  by = c("speciesName" = "acceptedName")) %>%
  left_join(manual_add, by = "label", suffix = c("", ".manual")) %>%
  mutate(
    type            = coalesce(type.manual, type),
    branchiness     = coalesce(branchiness.manual, branchiness),
    flowers         = coalesce(flowers.manual, flowers),
    avg_size        = coalesce(avg_size.manual, avg_size),
    density         = coalesce(density.manual, density),
    elevation_bin = coalesce(as.character(elevation_bin.manual), as.character(elevation_bin))
  ) %>%
  select(label, speciesName, type, branchiness, flowers, density, avg_size, elevation_bin) %>%
  arrange(label)

trait_table <- trait_table %>%
  mutate(
    speciesName = case_when(
      grepl("_sp__", label) ~ sub("^([^_]+)_sp__([^_]+).*", "\\1 sp \\2", label),
      grepl("_cf_", label) ~ sub("^([^_]+)_cf_([^_]+).*", "\\1 cf \\2", label),
      TRUE ~ sub("^([^_]+)_([^_]+).*", "\\1 \\2", label)
    )
  )

trait_table <- trait_table %>%
  mutate(
    flowers = round(flowers, 2),
    avg_size = round(avg_size, 2),
    density = round(density, 2)
  )

trait_table <- trait_table %>%
  rename(elevation = elevation_bin) %>%
  select(-elevation_state)     

write.csv(trait_table, "data/trait_table_all_species.csv", row.names = FALSE, quote = FALSE)
