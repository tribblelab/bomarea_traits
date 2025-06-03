setwd("~/Desktop/bomarea_traits/")

library(dplyr)
library(tidyr)
library(vegan)
library(readxl)
library(ape)
library(tidytree)
library(car)

traits <- read_xlsx("data_prep/bomarea_traits.xlsx", sheet = 1, na = "N/A")

#function
get_most_frequent_discrete_value <- function(vec) {
  tbl <- table(vec)
  if (dim(tbl) > 0) {
    value <- names(tbl)[which.max(tbl)]
    return(value)
  } else {
    return (NA)
  }
  }

### Name Check
sameName <- function(df) {
  if (is.na(df[2])) {
          return(NA)
  } else{
          if (df[1] != df[2]) {
                  return(1)
          } else {
                  return(NA)
          }   
  }
}

subsetName <- traits[, c("speciesName", "acceptedName")]
subsetName$ifSame <- apply(subsetName, 1, sameName)
subsetName <- na.omit(subsetName)
#write.csv(subsetName, file = "NameCheck")

### Name Check End

traits %>% 
  group_by(acceptedName) %>% 
  summarise(numBranchP = mean(numBranchP, na.rm = T), 
            numBracts = mean(numBracts, na.rm = T),
            numBranch1 = mean(numBranch1, na.rm = T), 
            numBranch2 = mean(numBranch2, na.rm = T),
            numBranch3 = mean(numBranch3, na.rm = T),
            numBranch4 = mean(numBranch4, na.rm = T),
            numBranch5 = mean(numBranch5, na.rm = T),
            Bracteoles = get_most_frequent_discrete_value(Bracteoles),
            ifFruiting = get_most_frequent_discrete_value(ifFruiting),
            matFruit = get_most_frequent_discrete_value(matFruit),
            ifFlowering = get_most_frequent_discrete_value(ifFlowering),
            matFlower = get_most_frequent_discrete_value(matFlower),
            allFlowersMat = get_most_frequent_discrete_value(allFlowersMat),
            nectarGuides = get_most_frequent_discrete_value(nectarGuides),
            colorTepalP = get_most_frequent_discrete_value(colorTepalP),
            colorTepalS = get_most_frequent_discrete_value(colorTepalS),
            ifTepalLengthMatch = get_most_frequent_discrete_value(ifTepalLengthMatch),
            ETepalLength = get_most_frequent_discrete_value(ETepalLength),
            ifExcerted = get_most_frequent_discrete_value(ifExcerted)) -> traits_by_species

traits_by_species_averaged <- cbind(traits_by_species[, c("acceptedName","numBranchP","numBracts","numBranch1")],
                                    traits_by_species[, c("Bracteoles","ifFruiting","matFruit","ifFlowering",       
                                                         "matFlower","allFlowersMat","nectarGuides",      
                                                         "colorTepalP","colorTepalS","ifTepalLengthMatch",
                                                         "ETepalLength","ifExcerted")]) #binding columns

#cleaning (inf to NA, Nan to NA, empty to NA, then converting those into NA instead of "NA")
traits_by_species_averaged[sapply(traits_by_species_averaged, is.infinite)] <- NA
traits_by_species_averaged[sapply(traits_by_species_averaged, is.nan)] <- NA
traits_by_species_averaged <- as.data.frame(apply(traits_by_species_averaged, 
                                                 2, 
                                                 car::recode, 
                                                 recodes = "'NA' = NA"))
traits_by_species_averaged <- as.data.frame(apply(traits_by_species_averaged, 
                                                  2, 
                                                  car::recode, 
                                                  recodes = "'' = NA"))
#turning these into numerical data
traitdata <- traits_by_species_averaged
traitdata$numBranchP <- as.numeric(traitdata$numBranchP)
traitdata$numBracts <- as.numeric(traitdata$numBracts)
traitdata$numBranch1 <- as.numeric(traitdata$numBranch1)


##############################
##### Inflorescence Type #####
##############################


#looking and 2 & 3 columns in dataset (branching and bracts)
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

traitdata %>%
  mutate(umbellike = numBranch1==0) %>% #adds new column TRUE if numBranch1 is 0, FALSE if other
  mutate(bracteoles = Bracteoles=="Y") %>% #adds new column TRUE if Bracteoles == "Y", FALSE if other
  mutate(acceptedName = gsub(" ", "_", acceptedName)) %>% #random but gets rid of "_" in accepted name
  select(acceptedName, umbellike, bracteoles) -> traitdatasubset #only keeps acceptedName, umbellike, bracteoles columns

traitdatasubset$type = apply(traitdatasubset, 1, typeset) #applies typeset function to data subset

###infl trait isolation (size + tepal traits)

traits %>%
  group_by(acceptedName) %>% #same as above for size and tepal traits
  summarise(maxBranchNo = max(numBranchP, na.rm = T),
            maxBranchLength = max(lengthTotal1, lengthTotal2, lengthTotal3, lengthTotal4, lengthTotal5, na.rm = T),
            degreeBranch = max(numBranch1, numBranch2, numBranch3, numBranch4, numBranch5, na.rm = T),
            colorTepalP = get_most_frequent_discrete_value(colorTepalP),
            colorTepalS = get_most_frequent_discrete_value(colorTepalS),
            ifTepalLengthMatch = get_most_frequent_discrete_value(ifTepalLengthMatch),
            ifExcerted = get_most_frequent_discrete_value(ifExcerted)) -> inflSelect #group and summarize
inflSelect %>%
  count(colorTepalP, sort = TRUE) #counts unique occurences and sorts
inflSelect %>%
  count(colorTepalS, sort = TRUE) #counts unique occurences and sorts


#taking the max num of flowers per species
traits$maxFlowers <- pmax(traits$numFlowers_1, traits$numFlowers_2, traits$numFlowers_3, na.rm = TRUE)

####################
##### Sparsity #####
####################

sparsity_df <- traits %>%
        mutate(numBranchesMeasured = rowSums(!is.na(select(., starts_with("lengthBranch"))))) %>%
        group_by(acceptedName) %>%
        slice_max(maxFlowers, with_ties = FALSE) %>%
        mutate(totalBranchLength = lengthTotal1 + 
                       coalesce(lengthSeg1_1, 0) + coalesce(lengthBranch1_1, 0) +
                       coalesce(lengthSeg1_2, 0) + coalesce(lengthBranch1_2, 0) +
                       coalesce(lengthSeg1_3, 0) + coalesce(lengthBranch1_3, 0) +
                       coalesce(lengthSeg1_4, 0) + coalesce(lengthBranch1_4, 0) +
                       coalesce(lengthSeg1_5, 0),
               sparsity = ifelse(numBranchesMeasured == 0, 0, ((totalBranchLength / numBranchesMeasured) * numBranchP) / 10)
        ) %>%
        select(gbfID, acceptedName, maxFlowers, numBranchP, numBranchesMeasured, totalBranchLength, sparsity) %>%
        ungroup()


#######################
##### Branchiness #####
#######################

branchiness_df <- traits %>%
        mutate(
                totNumBranches = rowSums(!is.na(select(., numBranch1, numBranch2, 
                                                        numBranch3, numBranch4, numBranch5)), na.rm = TRUE),
                branchiness = totNumBranches / numBranchP
        ) %>%
                select(gbfID, acceptedName, totNumBranches, numBranchP, branchiness) %>%
                ungroup()

        #averages of branchiness across species
        branchiness_max <- branchiness_df %>%
                group_by(acceptedName) %>%
                summarize(
                        branchiness = max(totNumBranches, na.rm = TRUE)
                ) %>%
                ungroup()

################
##### Size #####
################
        
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

######################
##### Tree Edits #####
######################
        
#match names to tree and drop some tips
tree <- read.tree("data/bom_only_MAP.tre")
tips_to_drop <- tree$tip.label[grep("caudata|herbertiana|glaucescens|Bomarea_edulis_Brazil_Campbell8900|Bomarea_edulis_Venezuela_Bunting4817|
                                        Bomarea_sp__catanatasoya_Peru_Graham12611|Bomarea_foliosa_Ecuador_Zak2268|Bomarea_lehmannii_AlzateS_N_
                                        |Bomarea_straminea_Alzate3300|Bomarea_pauciflora_AlzateS_N_|Bomarea_distichophylla_Peru_Rojas2689|
                                        Bomarea_superba_CultivatedinCA_SFBG20120052_A|Bomarea_foliosa_Ecuador_Zak2268|Bomarea_straminea_Alzate3300|
                                        Bomarea_lehmannii_AlzateS_N_|Bomarea_distichophylla_Peru_Rojas2689 ", 
                                        tree$tip.label)]
tree_edited <- ape::drop.tip(tree, tips_to_drop)
write.tree(tree_edited, file = "data/tree_edited.tre")
tree_df <- as_tibble(tree_edited)

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

#combine species names
tree_df$speciesName <- unlist(lapply(tree_df$label, get_gen_sp))
tree_df <- left_join(tree_df, traitdatasubset, by = c("speciesName" = "acceptedName"))
typedat <- data.frame(label = tree_df$label, type = tree_df$type)
typedat <- typedat[is.na(typedat$label) == FALSE, ]
typedat$type <- replace_na(as.character(as.integer(typedat$type)), "?")

#manually added these infl types to nexus file
manual_add <- data.frame(
  label = c("Bomarea_parvifolia_Peru_Stein2019", "Bomarea_tribachiata_AlzateS_N_", "Bomarea_angustipetala_Alzate5116",
            "Bomarea_lehmannii_AlzateS_N_", "Bomarea_straminea_Alzate3300", "Bomarea_chimborazensis_Ecuador_Aedo13023",
            "Bomarea_trimorphophylla_Alzate3158", "Bomarea_hartwegii_Alzate3157", "Bomarea_alstroemeriodes_Peru_Dillon1747",
            "Bomarea_euryphylla_Ecuador_Vargas2930", "Bomarea_foliosa_Ecuador_Zak2268", "Bomarea_bredemeyerana_Venezuela_Liesner7935",
            "Bomarea_acuminata_CR_Bonifacino6051", "Bomarea_bredemeyerana_Colombia_Tribble06", "Bomarea_killipii_Peru_Vasquez33143"),

  type = c(2, 2, 2, 2, 2, 2, 2, 0, 0, 1, 1, 0, 1, 0, 2)
)


#merge into existing df
typedat <- typedat %>%
  left_join(manual_add, by = "label") %>%
  mutate(type = ifelse(!is.na(type.y), type.y, type.x)) %>%
  select(label, type)

#turn all bracteoles into simple
typedat$type <- ifelse(typedat$type == "1", "0", typedat$type)

#make a matrix for type
typemat <- matrix(typedat$type, ncol =1)
rownames(typemat) <- typedat$label
colnames(typemat) <- "type"

#matrix for branchiness
branchiness_mat <- as.matrix(branchiness_max$branchiness)
rownames(branchiness_mat) <- branchiness_max$acceptedName
colnames(branchiness_mat) <- "branchiness"


#write to nexus 
write.nexus.data(typemat, file = "data/type_bi.nexus", format = "standard", missing = "?")
write.nexus.data(branchiness_mat, file = "data/branch.nexus", format = "standard", missing = "?")
