library(ape)
library(tidyverse)
setwd("~/Desktop/bomarea_traits")

source("scripts/latent_liability/print_latent_liability_functions.R")


### turn traits into df

# type
type_lines <- readLines("data/type.nexus")
start_type <- grep("MATRIX", type_lines, ignore.case = TRUE) + 1
end_type <- grep(";", type_lines[start_type:length(type_lines)], fixed = TRUE)[1] + start_type - 2

type_mat <- type_lines[start_type:end_type]
type_mat <- type_mat[!grepl("^\\s*$|\\[|#", type_mat)]

type_df <- read.table(text = type_mat, header = FALSE, stringsAsFactors = FALSE)
names(type_df) <- c("species", "type")

## continuous

# flowers
flowers_lines <- readLines("data/flowers.nexus")
start_flowers <- grep("MATRIX", flowers_lines, ignore.case = TRUE) + 1
end_flowers <- grep(";", flowers_lines[start_flowers:length(flowers_lines)], fixed = TRUE)[1] + start_flowers - 2

flowers_mat <- flowers_lines[start_flowers:end_flowers]
flowers_mat <- flowers_mat[!grepl("^\\s*$|\\[|#", flowers_mat)]

flowers_df <- read.table(text = flowers_mat, header = FALSE, stringsAsFactors = FALSE)
names(flowers_df) <- c("species", "flowers")

# size
size_lines <- readLines("data/size.nexus")
start_size <- grep("MATRIX", size_lines, ignore.case = TRUE) + 1
end_size <- grep(";", size_lines[start_size:length(size_lines)], fixed = TRUE)[1] + start_size - 2

size_mat <- size_lines[start_size:end_size]
size_mat <- size_mat[!grepl("^\\s*$|\\[|#", size_mat)]

size_df <- read.table(text = size_mat, header = FALSE, stringsAsFactors = FALSE)
names(size_df) <- c("species", "size")

# density
density_lines <- readLines("data/density.nexus")
start_density <- grep("MATRIX", density_lines, ignore.case = TRUE) + 1
end_density <- grep(";", density_lines[start_density:length(density_lines)], fixed = TRUE)[1] + start_density - 2

density_mat <- density_lines[start_density:end_density]
density_mat <- density_mat[!grepl("^\\s*$|\\[|#", density_mat)]

density_df <- read.table(text = density_mat, header = FALSE, stringsAsFactors = FALSE)
names(density_df) <- c("species", "density")

### Make df for .xml file
## Load in data and tree (D,C,C)

# 3 variables
bom.data <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE),
                   list(type_df, size_df, density_df))
bom.tree <- read.tree("data/tree_edited.tre")

bom.data.for.xml <- bom.data[, -1]
names(bom.data.for.xml) <- c("d", "c", "c")
row.names(bom.data.for.xml) <- bom.data$species

bom.data.for.xml[, 1] <- gsub("\\?", NA, bom.data.for.xml[, 1])
bom.data.for.xml[, 1] <- as.numeric(bom.data.for.xml[, 1])

# bom.data.for.xml <- bom.data.for.xml[, c(1, 3, 4)]
# names(bom.data.for.xml) <- c("d", "c", "c")

bom.data.for.xml[, 2:3] <- lapply(bom.data.for.xml[, 2:3], as.numeric)
bom.data.for.xml[, 2:3] <- scale(bom.data.for.xml[, 2:3])

# 2 variables
bom.data <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE),
                   list(type_df, density_df))
bom.tree <- read.tree("data/tree_edited.tre")

bom.data.for.xml <- bom.data[, -1]
names(bom.data.for.xml) <- c("d", "c")
row.names(bom.data.for.xml) <- bom.data$species

bom.data.for.xml[, 1] <- gsub("\\?", NA, bom.data.for.xml[, 1])
bom.data.for.xml[, 1] <- as.numeric(bom.data.for.xml[, 1])

# bom.data.for.xml <- bom.data.for.xml[, c(1, 3)]
# names(bom.data.for.xml) <- c("d", "c")

bom.data.for.xml[, 2] <- scale(bom.data.for.xml[, 2])

## Write the xml
printLatentLiability(file="scripts/latent_liability/size_flowers_log/latent_liability_size_flowers_log.xml",
                     latent.liability.info=bom.data.for.xml,
                     tree=bom.tree,
                     log.name="bomarea_latent_liability_size_flowers_log",
                     ngen="100000000", log.every="10000",
                     walk.or.scale="walk",
                     is.multistate = 1,
                     name.for.traits="bomareaTraits",
                     jitter=0, precision=0.05, wishart=0.1)
