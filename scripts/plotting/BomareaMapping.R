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
