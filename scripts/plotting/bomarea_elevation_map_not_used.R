#get Bomarea data using BIEN network, get altitude data for each individual, and make boxplots of species 
library("jsonlite")
library("rgbif")
library(maptools)
library(sp)
library(scales)
library(tidyverse)
setwd("~/Documents/Berkeley_IB/research_projects/Bomarea/bomarea_dist_alt")

# load in world shapefile
data("wrld_simpl")

##############
read.csv("bomarea_alt_all.csv",stringsAsFactors=F ) %>%
  na.omit() -> bomarea
  
#############

coordinates(bomarea) <- ~longitude + latitude
proj4string(bomarea) <- wrld_simpl@proj4string

plot(bomarea, pch = 16)

rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "aliceblue")

plot(wrld_simpl, col="ivory", add = T)

plot(bomarea, pch=16, col = alpha("black", 0.4), add = T)





species<-unique(bomarea$scrubbed_species_binomial)
length(unique(bomarea$scrubbed_species_binomial))
species<-data.frame(species, c(1:124))
colnames(species)<-c("scrubbed_species_binomial","color")
bomarea<-merge(bomarea, species, by="scrubbed_species_binomial")

bomarea<-bomarea[!is.na(bomarea$longitude),]


##### Now add elevation


#### create subset of BIEN network data for the species in Antioquia
## import list of species in antioquia
bom_ant<-read.csv("BomareaInAntioquia.csv")
spAnt<-bom_ant$Scientific.Name

bomareaCol <- bomarea[which(bomarea[,1] %in% spAnt),]


##### the max requests per day is 2500
##### There is no rate max. 

apikey="AIzaSyCHplAQGySUgGcHr6cOOFmUKCeyo3LQgGo"

#bomareaCol<-head(bomareaCol)



output<-data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(output)<-c("latitude","longitude","elevation","species")
output$species<-as.character(output$species)
for (i in 1:length(bomareaCol$latitude))
{
coord<-data.frame(latitude=bomareaCol$latitude[i], longitude=bomareaCol$longitude[i])
names(coord)<-c("decimalLatitude","decimalLongitude")
A<-elevation(input=coord, key=apikey)  
B<-cbind(A,bomareaCol$scrubbed_species_binomial[i])
colnames(B)<-c("latitude","longitude","elevation","species")
output<-rbind(output,B)
  
}

write.csv(output,file="bomaraColalt.csv")

bomareaCollist<-split(output, output$species)
spAnt<-names(bomareaCollist)


means<-NA
for (i in 1:length(bomareaCollist)) 
{
  means[i]<-mean(bomareaCollist[[i]]$elevation, na.rm=T) 
}
means<- data.frame(means,spAnt)
names(means)<-c("meanalt","scrubbed_species_binomial")
colnames(output)<-c("latitude", "longitude", "elevation","scrubbed_species_binomial")
bomareaCol2<-merge(output,means,by="scrubbed_species_binomial")

bomareaCol3<-bomareaCol2[order(bomareaCol2$meanalt),]

#boxplot(bomarea3$alt$elevation, main="Bomarea Altitude Distribution", ylab="meters")

p <- ggplot(bomareaCol3, aes(factor(meanalt), elevation))

pdf("Bomarea_Antioquia_Altitudes.pdf")
p + geom_boxplot()  + 
  geom_boxplot() +  
  labs(x="Species", y="Altitude (m)",title = "Altitude Distribution in Bomarea Species in Antioquia") +
  scale_x_discrete(labels=unique(bomareaCol3$scrubbed_species_binomial)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, color="black"))
dev.off()

write.csv(bomareaCol3,file="FinalAntBomareaAltData.csv")


###### for all bomareas, save data 

output<-data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(output)<-c("latitude","longitude","elevation","species")
output$species<-as.character(output$species)
for (i in 1:length(bomarea$latitude))
{
  coord<-data.frame(latitude=bomarea$latitude[i], longitude=bomarea$longitude[i])
  names(coord)<-c("decimalLatitude","decimalLongitude")
  A<-elevation(input=coord, key=apikey)  
  B<-cbind(A,bomarea$scrubbed_species_binomial[i])
  colnames(B)<-c("latitude","longitude","elevation","species")
  output<-rbind(output,B)
  
}



write.csv(output,file="bomarea_alt_all.csv")
bomarea_alt_all<-read.csv("bomarea_alt_all.csv")

bomarea_alt_all<-bomarea_alt_all[-1,]
bomarea_alt_all<-bomarea_alt_all[,-1]

bomarea_spl<-split(bomarea_alt_all, bomarea_alt_all$species)
names(bomarea_spl)<-na.omit(unique(bomarea_alt_all$species))

##get means 


means<-NA
for (i in 1:length(bomarea_spl)) 
{
  means[i]<-mean(bomarea_spl[[i]]$elevation, na.rm=T) 
}
means<- data.frame(means,names(bomarea_spl))

bomarea_means<-means[order(means$means),]
bomarea_means$index=c(1:length(bomarea_means$means))

plot(bomarea_means$means,bomarea_means$index, pch=19, 
     xlab="Mean Altitudes of Bomarea spp. by Species (m)", yaxt='n', ylab=NA)
## find model to fit to this line 

plot((bomarea_means$index)^(1/3), bomarea_means$means, pch=19, 
     xlab="Mean Altitudes of Bomarea spp. by Species (m)", yaxt='n', ylab=NA)
a<-lm(bomarea_means$means~bomarea_means$index^(1/3))


# ── Bomarea downloaded elevation plot from coordinates only ───────────────────

library(ape)
library(dplyr)
library(ggplot2)
library(sf)
library(elevatr)

setwd("~/Desktop/bomarea_traits/")

# ── 1. Species list from NEXUS ────────────────────────────────────────────────

nex <- read.nexus.data("data/type.nexus")
taxa_names <- names(nex)

species_raw <- sapply(strsplit(taxa_names, "_"),
                      function(x) paste(x[1], x[2]))

species_clean <- unique(
  species_raw[
    !grepl("^\\S+ (cf|sp)$",
            species_raw,
            ignore.case = TRUE)
  ]
)

cat("Species retained:", length(species_clean), "\n")

# ── 2. Load GBIF occurrence data ──────────────────────────────────────────────

occ <- read.delim(
  "data_prep/occurence_data.csv",
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)

# ── 3. Filter coordinates only ────────────────────────────────────────────────

occ_filtered <- occ %>%
  filter(
    species %in% species_clean,
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),

    decimalLongitude >= -125,
    decimalLongitude <= -30,

    decimalLatitude >= -55,
    decimalLatitude <= 45,

    # remove California cultivated-region records
    !(decimalLongitude >= -125 &
        decimalLongitude <= -113 &
        decimalLatitude >= 32 &
        decimalLatitude <= 42)
  )

cat("Occurrences retained:", nrow(occ_filtered), "\n")

# ── 4. Convert coordinates to sf object ───────────────────────────────────────

occ_sf <- st_as_sf(
  occ_filtered,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE
)

# ── 5. DOWNLOAD elevations from coordinates only ──────────────────────────────

occ_elev <- get_elev_point(
  locations = occ_sf,
  src = "aws",
  z = 7,
  prj = "EPSG:4326",
  overwrite = TRUE
)

cat("Elevations recovered:",
    sum(!is.na(occ_elev$elevation)), "\n")

cat("Missing elevations:",
    sum(is.na(occ_elev$elevation)), "\n")

# ── 6. Clean elevation dataframe ──────────────────────────────────────────────

occ_elev_df <- occ_elev %>%
  st_drop_geometry() %>%
  rename(elevation_m = elevation) %>%
  filter(!is.na(elevation_m))

# ── 7. Calculate species mean elevations ──────────────────────────────────────

species_means <- occ_elev_df %>%
  group_by(species) %>%
  summarise(
    mean_elev = mean(elevation_m, na.rm = TRUE),
    min_elev = min(elevation_m, na.rm = TRUE),
    max_elev = max(elevation_m, na.rm = TRUE),
    n_records = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_elev)

species_means$index <- 1:nrow(species_means)

# Remove negative elevations first
occ_elev_df <- occ_elev_df %>%
  filter(elevation_m >= 0)

# Recalculate species means after removing negatives
species_means <- occ_elev_df %>%
  group_by(species) %>%
  summarise(
    mean_elev = mean(elevation_m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_elev) %>%
  mutate(index = row_number())

# Read type.nexus
type_nex <- read.nexus.data("data/type.nexus")

type_df <- data.frame(
  label = names(type_nex),
  type_state = as.character(unlist(type_nex)),
  stringsAsFactors = FALSE
)

type_df$species <- sapply(
  strsplit(type_df$label, "_"),
  function(x) paste(x[1], x[2])
)

type_df <- type_df %>%
  filter(!grepl("^\\S+ (cf|sp)$", species, ignore.case = TRUE))

type_df$type <- ifelse(
  type_df$type_state == "0",
  "Umbel",
  "Compound"
)

type_df <- type_df %>%
  distinct(species, type)

# Join type ONCE
species_means <- species_means %>%
  left_join(type_df, by = "species") %>%
  arrange(mean_elev)

# Set species order by increasing mean elevation
species_order <- species_means$species

occ_elev_df <- occ_elev_df %>%
  left_join(type_df, by = "species") %>%
  mutate(
    species = factor(species, levels = species_order)
  )

species_means <- species_means %>%
  mutate(
    species = factor(species, levels = species_order)
  )

# Order species by mean elevation
species_order <- occ_elev_df %>%
  group_by(species) %>%
  summarise(
    mean_elev = mean(elevation_m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_elev) %>%
  pull(species)

# Apply ordering
occ_elev_df$species <- factor(
  occ_elev_df$species,
  levels = species_order
)

# Boxplot ordered by increasing elevation and colored by type
p_box <- ggplot(occ_elev_df,
                aes(x = species,
                    y = elevation_m,
                    fill = type)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "Species ordered by mean elevation",
    y = "Downloaded elevation (m)",
    fill = "Inflorescence type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 5,
      colour = "black"
    )
  )

ggsave(
  "figures/bomarea_elevation_boxplot_by_type.pdf",
  plot = p_box,
  width = 14,
  height = 6,
  device = "pdf"
)

ggsave(
  "figures/bomarea_elevation_boxplot_by_type.png",
  plot = p_box,
  width = 14,
  height = 6,
  dpi = 300
)
