library(leaflet)
library(readxl)
library(writexl)
library(dplyr)
library(RColorBrewer)
library(htmlwidgets)
library(ape)

setwd("~/Desktop/bomarea_traits/data/")

#######################################
##data cleaning

#adding inflorescence type to spreadsheet
nex_file <- read.nexus.data("type_dropped_tips.nexus")
nex_df <- data.frame(
  Species = names(nex_file),
  Inflorescence_type = as.numeric(unlist((nex_file)))
)
nex_df <- nex_df[!is.na(nex_df$Inflorescence_type), ]

#merge the data with the spreadsheet (only do once)
xl <- read_excel("geo_data.xlsx")
geo_data <- xl %>%
  left_join(nex_df, by = "Species")

#write_xlsx(geo_data, "geo_data.xlsx")
convert_to_decimal <- function(coord) {
  #remove brackets if they exist
  coord <- gsub("[\\[\\]]", "", coord)

  #match degrees and minutes
  matches <- regmatches(coord, gregexpr("\\d+", coord))
  d <- as.numeric(matches[[1]][1]) # Degrees
  m <- as.numeric(matches[[1]][2]) # Minutes

  #decimal degrees
  decimal <- d + m / 60

  #check for direction and reassign
  if (grepl("[SW]", coord)) {
    decimal <- -decimal
  }

  return(decimal)
}

#apply to lat and long
geo_data$Latitude <- sapply(geo_data$Latitude, convert_to_decimal)
geo_data$Longitude <- sapply(geo_data$Longitude, convert_to_decimal)

#######################################
##make map

#load geo_data
geo_data <- read_xlsx("geo_data.xlsx")
geo_data <- geo_data[!is.na(geo_data$Inflorescence_type), ]

#using RColorBrewer
num_colors <- length(unique(geo_data$Inflorescence_type)) 
colors <- brewer.pal(num_colors, "Set1")

# Assign colors to the Inflorescence_type
geo_data$Type_color <- colors[as.factor(geo_data$Inflorescence_type)]

#create map
map <- leaflet(geo_data) %>%
  addTiles() %>%
  addCircleMarkers(
    lat = ~Latitude,
    lng = ~Longitude,
    popup = ~paste("Species: ", Species,
                  "<br>Inflorescence Type: ", Inflorescence_type),
    radius = 5,
    color = ~Type_color,
    stroke = FALSE,
    fillOpacity = 0.7
  )
map

saveWidget(map, "infl_type_map.html")
