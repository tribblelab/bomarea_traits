library(ape)
library(tidyverse)
source("beast/print_latent_liability_xml(1).R")

### Turn traits into df

# Type
type_lines <- readLines("bomarea_traits/data/binary_type.nexus")
start_type <- grep("MATRIX", type_lines, ignore.case = TRUE) + 1
end_type <- grep(";", type_lines[start_type:length(type_lines)], fixed = TRUE)[1] + start_type - 2

type_mat <- type_lines[start_type:end_type]
type_mat <- type_mat[!grepl("^\\s*$|\\[|#", type_mat)]

type_df <- read.table(text = type_mat, header = FALSE, stringsAsFactors = FALSE)
names(type_df) <- c("species", "type")
type_df$Species <- trimws(type_df$species)

## Continuous

# Size
size_lines <- readLines("bomarea_traits/data/size.nexus")
start_size <- grep("MATRIX", size_lines, ignore.case = TRUE) + 1
end_size <- grep(";", size_lines[start_size:length(size_lines)], fixed = TRUE)[1] + start_size - 2

size_mat <- size_lines[start_size:end_size]
size_mat <- size_mat[!grepl("^\\s*$|\\[|#", size_mat)]

size_df <- read.table(text = size_mat, header = FALSE, stringsAsFactors = FALSE)
names(size_df) <- c("species", "size")

# Sparsity
sparsity_lines <- readLines("bomarea_traits/data/sparsity.nexus")
start_sparsity <- grep("MATRIX", sparsity_lines, ignore.case = TRUE) + 1
end_sparsity <- grep(";", sparsity_lines[start_sparsity:length(sparsity_lines)], fixed = TRUE)[1] + start_sparsity - 2

sparsity_mat <- sparsity_lines[start_sparsity:end_sparsity]
sparsity_mat <- sparsity_mat[!grepl("^\\s*$|\\[|#", sparsity_mat)]

sparsity_df <- read.table(text = sparsity_mat, header = FALSE, stringsAsFactors = FALSE)
names(sparsity_df) <- c("species", "sparsity")

### Make df for .xml file
## Load in data and tree (D,C,C)
bom.data <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE),
                   list(type_df, size_df, sparsity_df))
bom.tree <- read.tree("bomarea_traits/data/tree_edited.tre")

bom.data.for.xml <- bom.data[, -1]
names(bom.data.for.xml) <- c("d", "c", "c")
row.names(bom.data.for.xml) <- bom.data$species

bom.data.for.xml[, 1] <- gsub("\\?", NA, bom.data.for.xml[, 1])
bom.data.for.xml[, 1] <- as.numeric(bom.data.for.xml[, 1])

bom.data.for.xml <- bom.data.for.xml[, c(1, 3, 4)]
names(bom.data.for.xml) <- c("d", "c", "c")

bom.data.for.xml[, 2:3] <- scale(bom.data.for.xml[, 2:3])

## Write the xml
printLatentLiability(file="bomarea_traits/scripts/latent_liability/bomarea_latent_liability.xml",
                     latent.liability.info=bom.data.for.xml,
                     tree=bom.tree,
                     log.name="bomarea_latent_liability",
                     ngen="100000000", log.every="10000",
                     walk.or.scale="walk",
                     name.for.traits="bomareaTraits",
                     jitter=0.2, precision=0.05, wishart=0.1,
                     is.multistate=1)

