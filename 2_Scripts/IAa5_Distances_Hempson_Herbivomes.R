######################################################################!
#                                                                    #
#      IAa5: SELECT NON-ANALOGUE CELLS                               #  
#                                                                    #
######################################################################!

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # METADATA
    # This script loads the earlier saved Hempson distance matrix and determines how dissimilar the different Hempson herbivomes are.
    #
    # Output: a table indicating giving the median Squared Chord Distances between and within the cells of different Hempson Herbivomes.
    # time: needs 1 min to run.
    # Status: Ready (2021/01/05)
    # Note: This scrips used large matrices (2-10 Gb). This is heavier than R can work with. Therefore we use the Bigmemory package.
    #       This package requires loading the csv file distance matrix. This takes long the first time, but should be quicker afterwards.
    #
    #         I ran this script on 2025/02/03
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 1. PREPARING THE DATA ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Clear Memory ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

rm(list = ls())  # clean memory (= remove)
graphics.off()  # close graphic windows

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Set working directory ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

getwd()


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Set up libraries ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

library(dplyr) # for data cleaning

# To calculate the distances
#–––––––––––––––––––––––––––––––––––––––––––––––––
library(analogue) # to do the cut-off ROC analysis
library(pheatmap) # to make heathmaps

# to work with large datasets
#–––––––––––––––––––––––––––––––––––––––––––––––––
library(bigmemory) # To work with large datasets within R
library(biganalytics)
library(bigalgebra)

# for GIS stuff
#–––––––––––––––––––––––––––––––––––––––––––––––––
library(sp) # GIS
library(rgeos) # GIS
library(rgdal) # GIS
library(raster)
library(ggplot2)
library(ggfortify)
library("rworldmap") # visualisation data

# For plotting
#–––––––––––––––––––––––––––––––––––––––––––––––––
library(ggpubr)
library(viridis)
library(gridExtra) #grid arrange
library(RColorBrewer) # colour scheme
library(scales) # to scale colours in ggplot2


## – Start Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
start_time <- Sys.time() # Save Time


## – Directories  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
if(!dir.exists(paste0("./3_Output/IAa5_Calculate_Distances_Herbivomes/"))){
  dir.create(paste0("./3_Output/IAa5_Calculate_Distances_Herbivomes/"))}



#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Load Data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
functional.community.matrix.relative = read.csv("./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_Hempson_relative_abundances.csv")
functional.community.matrix.relative = functional.community.matrix.relative[,!(colnames(functional.community.matrix.relative) %in% c("X","CellNr"))] # remove two redundant variables

# Spatial
#–––––––––––––––––––––––––––––––––––––––––––––––––
raster.example = raster("~/Datasets/PHYLACINE_V1_2_1/Ranges/Present_natural/Stegodon_orientalis.tif")
biogeo.realms = readOGR("/Users/au572919/Datasets/Biogeographical_Realms/udvardy/udvardy_py.shp")
herbivomes = readOGR("./1_files/Herbivome_Polygons.shp")

# Attach (= load) the big distance matrix
#–––––––––––––––––––––––––––––––––––––––––––––––––
big.distance.matrix <- attach.big.matrix("./3_Output/IAa4_Calculate_Cell_Distances/IAa4_big_distance_matrix_Hempson.desc")

# Check data
#–––––––––––––––––––––––––––––––––––––––––––––––––
big.distance.matrix[1:10,1:10]
dim(big.distance.matrix)
# range(big.distance.matrix[,]) #range = 0-200

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Load Functions ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

source("./2_scripts/Scripts_2020_12_10/Functions/I_F1_Realm_Distances.R")

# Define function 
#–––––––––––––––––––––––––––––––––––––––––––––––––
mean.without.na = function(x) {mean(na.omit(x))}
max.without.na = function(x) {max(na.omit(x))}
median.without.na = function(x) {median(na.omit(x))}
min.without.na = function(x) {min(na.omit(x))}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# –  Prepare spatial data ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Make Countries Polygon
#–––––––––––––––––––––––––––––––––––––––––––––––––
data("countriesCoarse")
countries = (countriesCoarse)
countries@data$continent[which(countries@data$NAME == "French Guiana")] = "South America"
countries = spTransform(countries, CRSobj = crs(raster.example))
countries.poly = (countries[(countries@data$continent %in% c("Eurasia", "Africa", "Australia", "North America", "South America")),])

countries.raster = rasterize(countries.poly,raster.example, field = 1)
#plot(countries.raster)
#plot(countries.poly, add = T)

# Change the CRS 
#–––––––––––––––––––––––––––––––––––––––––––––––––
herbivomes.reprojected = spTransform(herbivomes, proj4string(raster.example))

# Change rasters values == 0 to NA's 
#–––––––––––––––––––––––––––––––––––––––––––––––––
raster.example[raster.example >= 0] <- NA


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Create data frames ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

df.master = functional.community.matrix.relative[,c("Longitude", "Latitude","Scenario", "Original.Cell", "Biogeo.realm", "Countries.raster", "Herbivores.Present")]
df.master = df.master[which(df.master$Herbivores.Present == 1),]

df.Hempson = df.master[df.master$Scenario == "Hempson",]

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 3. DIVIDE DISTANCE MATRIX ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# Divide cells between current and present-natural cells
Cells.Hempson = which(functional.community.matrix.relative[which(functional.community.matrix.relative$Herbivores.Present == 1),"Scenario"] == "Hempson")

      #~~~~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 
      #     Squared chord distance values can range from 0.0 to 2.0, with 0.0 indicating identical proportions of species within the samples being compared. 
      #     For the squared chord-distance metric, critical values have ranged from 0.05 in local applications to 0.40 for continental-scale paleoclimate reconstructions.
      #     values of 0.12 to 0.20 most commonly used in vegetational and paleoclimatic inference at intermediate scales.
      #     Several authors have used multiple critical values to indicate varying degrees of confidence in analog determination or to draw paleoenvironmental inferences at multiple ecological scales 
      #
      #     SOURCE:  Jackson, S. T., and J. W. Williams. 2004. MODERN ANALOGS IN QUATERNARY PALEOECOLOGY: Here Today, Gone Yesterday, Gone Tomorrow? Annual Review of Earth and Planetary Sciences 32:495–537.
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Subset Data (distance matrix) ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
nrow(functional.community.matrix.relative)

# How similar are current cells to each other?
#–––––––––––––––––––––––––––––––––––––––––––––––––
Hempson.vs.Hempson <- sub.big.matrix(big.distance.matrix,
                                     firstRow = 1,
                                     lastRow = length(Cells.Hempson),
                                     firstCol = 1,
                                     lastCol = length(Cells.Hempson))
dim(Hempson.vs.Hempson)

## Turn Diagonal into NA
#–––––––––––––––––––––––––––––––––––––––––––––––––
Hempson.vs.Hempson.with.NA= Hempson.vs.Hempson
d <- min(dim(Hempson.vs.Hempson.with.NA))
Hempson.vs.Hempson.with.NA[cbind(1:d, 1:d)] = NaN # Is the same as NA, but in C++

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 4. DISTANCE BETWEEN HEMPSOM HERBIVOMES ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

      #~~~~~ NOTE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #       
      #   This section calculates how large the distances are between the realms defined by Hempson.
      #   This value will be used as a threshold against which other changes can be compared.
      #   Note that we first used the current maps to calculate the herbivome differences, but the ranges of four species (see Hempson paper) have been changed to the present-natural. 
      #   This was done because they better resemble the historical ranges estimates by Hempson, which supposedly portray what herbivore assemblages were like 1000 years ago, i.e. before significant human impact.
      # 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create buffer (currently 2*50 km)
#–––––––––––––––––––––––––––––––––––––––––––––––––
plot(herbivomes.reprojected)
herbivomes.reprojected.buffer = gBuffer(herbivomes.reprojected, byid=T, id=NULL,
                                        width=-150*10^3,
                                        quadsegs=5, capStyle="ROUND",joinStyle="ROUND", mitreLimit=1.0)
plot(herbivomes.reprojected.buffer)

## Change data type from character to numeric (to avoid errors)
herbivomes.reprojected.buffer$Herbivome = as.numeric(herbivomes.reprojected.buffer$Herbivome)

# –  Transform the polygon into a raster
#–––––––––––––––––––––––––––––––––––––––––––––––––
herbivomes.reprojected.raster = rasterize(herbivomes.reprojected.buffer,
                                          raster.example,
                                          field = "Herbivome",
                                          fun='last',
                                          background=NA,
                                          mask=F,
                                          update=T,
                                          na.rm=TRUE)

#plot(herbivomes.reprojected.raster)

## Turn rastervalues into vector
herbivome.vector = getValues(herbivomes.reprojected.raster)

# Select the cells where herbivores occur
cells.with.herbivome.classification = which(!(is.na(herbivome.vector)))

# Identify which herbivore assemblages belong to each of the herbivomes
cells.with.herbivome.classification.Hempson = cells.with.herbivome.classification[which(cells.with.herbivome.classification %in% df.Hempson$Original.Cell)]

# Create new variable denoting the herbivome classification
df.Hempson$Herbivome.Hempson = NA

df.Hempson[which(df.Hempson$Original.Cell %in% cells.with.herbivome.classification.Hempson),"Herbivome.Hempson"] = herbivome.vector[cells.with.herbivome.classification.Hempson]


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Calculate the distance to other herbivomes ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––


## Rainforest
#–––––––––––––––––––––––––––––––––––––––––––––––––
african.rainforest.herbivome = which(df.Hempson$Herbivome.Hempson == 1)
not.african.rainforest.herbivome = which(df.Hempson$Herbivome.Hempson %in% c(2,3,4))

Median.1.between = median(apply(Hempson.vs.Hempson.with.NA[african.rainforest.herbivome, not.african.rainforest.herbivome],1,median.without.na))
Mean.1.between = mean(apply(Hempson.vs.Hempson.with.NA[african.rainforest.herbivome, not.african.rainforest.herbivome],1,mean.without.na))
Median.1.within = median(apply(Hempson.vs.Hempson.with.NA[african.rainforest.herbivome, african.rainforest.herbivome],1,median.without.na))


## Desert
#–––––––––––––––––––––––––––––––––––––––––––––––––

african.desert.herbivome = which(df.Hempson$Herbivome.Hempson == 2)
not.african.desert.herbivome = which(df.Hempson$Herbivome.Hempson %in% c(1,3,4))

Median.2.between = median(apply(Hempson.vs.Hempson.with.NA[african.desert.herbivome, not.african.desert.herbivome],1,median.without.na))
Mean.2.between = mean(apply(Hempson.vs.Hempson.with.NA[african.desert.herbivome, not.african.desert.herbivome],1,mean.without.na))
Median.2.within = median(apply(Hempson.vs.Hempson.with.NA[african.desert.herbivome, african.desert.herbivome],1,median.without.na))


## Savanna A (High VALS)
#–––––––––––––––––––––––––––––––––––––––––––––––––
african.savannahA.herbivome = which(df.Hempson$Herbivome.Hempson == 3)
not.african.savannahA.herbivome = which(df.Hempson$Herbivome.Hempson %in% c(1,2,4))

Median.3.between = median(apply(Hempson.vs.Hempson.with.NA[african.savannahA.herbivome, not.african.savannahA.herbivome],1,median.without.na))
Mean.3.between = mean(apply(Hempson.vs.Hempson.with.NA[african.savannahA.herbivome, not.african.savannahA.herbivome],1,mean.without.na))
Median.3.within = median(apply(Hempson.vs.Hempson.with.NA[african.savannahA.herbivome, african.savannahA.herbivome],1,median.without.na))



## Savanna B(Bulk Feeder)
#–––––––––––––––––––––––––––––––––––––––––––––––––
african.savannahB.herbivome = which(df.Hempson$Herbivome.Hempson == 4)
not.african.savannahB.herbivome = which(df.Hempson$Herbivome.Hempson %in% c(1,2,3))

Median.4.between = median(apply(Hempson.vs.Hempson.with.NA[african.savannahB.herbivome, not.african.savannahB.herbivome],1,median.without.na))
Mean.4.between = mean(apply(Hempson.vs.Hempson.with.NA[african.savannahB.herbivome, not.african.savannahB.herbivome],1,mean.without.na))
Median.4.within = median(apply(Hempson.vs.Hempson.with.NA[african.savannahB.herbivome, african.savannahB.herbivome],1,median.without.na))


## Savanna A & B
#–––––––––––––––––––––––––––––––––––––––––––––––––
african.savanna.herbivome = which(df.Hempson$Herbivome.Hempson %in% c(3,4))
not.african.savanna.herbivome  = which(df.Hempson$Herbivome.Hempson %in% c(1,2))

Median.3.4.between = median(apply(Hempson.vs.Hempson.with.NA[african.savanna.herbivome, not.african.savanna.herbivome],1,median.without.na))
Mean.3.4.between = mean(apply(Hempson.vs.Hempson.with.NA[african.savanna.herbivome, not.african.savanna.herbivome],1,mean.without.na))
Median.3.4.within = median(apply(Hempson.vs.Hempson.with.NA[african.savanna.herbivome, african.savanna.herbivome],1,median.without.na))

# distances between Hempson's herbivore communities
#–––––––––––––––––––––––––––––––––––––––––––––––––
round(c(Median.1.between, Median.2.between, Median.3.between, Median.4.between, Median.3.4.between),2)
round(c(Median.1.within, Median.2.within, Median.3.within, Median.4.within, Median.3.4.within),2)


## Create df
breaks = seq(min(Hempson.vs.Hempson.with.NA[african.rainforest.herbivome, not.african.rainforest.herbivome]), max(Hempson.vs.Hempson.with.NA[african.rainforest.herbivome, not.african.rainforest.herbivome]),by = 0.05)
col.1 <- rgb(0, 0, 225, max = 255, alpha = 100)
col.2 <- rgb(225, 0, 0, max = 255, alpha = 100)
ylim = c(0,20000)
#ylim = c(0,5000)


breaks = 50
similarities.df = as.data.frame(matrix(NA,5,2))
colnames(similarities.df) = c("between","within")
rownames(similarities.df) = c("Forest Duiker", "Arid Gazelle", "High VALS", "Bulk Feeder", "Savanna Combined")

png(filename = "./3_Output/IAa5_Calculate_Distances_Herbivomes/IAa5_Distances_herbivomes.png", width = 1200, height = 600)
par(mfrow=c(1,5))
hist(Hempson.vs.Hempson.with.NA[african.rainforest.herbivome, african.rainforest.herbivome], breaks = breaks, xlim = c(0,2), main = "Forest Duiker", col = col.2, ylim = ylim)
hist(Hempson.vs.Hempson.with.NA[african.rainforest.herbivome, not.african.rainforest.herbivome], breaks = breaks, xlim = c(0,2), col = col.1 , add = T)

hist(Hempson.vs.Hempson.with.NA[african.desert.herbivome, african.desert.herbivome], breaks = breaks, xlim = c(0,2), main = "Arid Gazelle", col = col.2, ylim = ylim)
hist(Hempson.vs.Hempson.with.NA[african.desert.herbivome, not.african.desert.herbivome], breaks = breaks, xlim = c(0,2), col = col.1, add = T)

hist(Hempson.vs.Hempson.with.NA[african.savannahA.herbivome, african.savannahA.herbivome], breaks = breaks, xlim = c(0,2), main = "High VALS", col = col.2, ylim = ylim)
hist(Hempson.vs.Hempson.with.NA[african.savannahA.herbivome, not.african.savannahA.herbivome], breaks = breaks, xlim = c(0,2), col = col.1, add = T)

hist(Hempson.vs.Hempson.with.NA[african.savannahB.herbivome, african.savannahB.herbivome], breaks = breaks, xlim = c(0,2), main = "Bulk Feeder", col = col.2, ylim = ylim)
hist(Hempson.vs.Hempson.with.NA[african.savannahB.herbivome, not.african.savannahB.herbivome], breaks = breaks, xlim = c(0,2), col = col.1, add = T)

hist(Hempson.vs.Hempson.with.NA[african.savanna.herbivome, african.savanna.herbivome], breaks = breaks, xlim = c(0,2), main = "Savanna Combined", col = col.2, ylim = ylim)
hist(Hempson.vs.Hempson.with.NA[african.savanna.herbivome, not.african.savanna.herbivome], breaks = breaks, xlim = c(0,2), col = col.1, add = T)
dev.off()




## Fill df
similarities.df[,"between"] = round(c(Median.1.between, Median.2.between, Median.3.between, Median.4.between, Median.3.4.between),2)
similarities.df[,"within"] = round(c(Median.1.within, Median.2.within, Median.3.within, Median.4.within, Median.3.4.within),2)

      #~~~~~ NOTE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #       
      #   WITHIN: the median distance between a cell and all cells within a herbivome.
      #   BETWEEN: the distance between a cell and all cells of different herbivomes
      #
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

similarities.df$ratio = round(similarities.df$between / similarities.df$within,2) # we want the ratios to be high

write.csv(similarities.df, "./3_Output/IAa5_Calculate_Distances_Herbivomes/IAa5a_Distance_Hempson_Herbivomes.csv")




## – End Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
end_time <- Sys.time()
end_time - start_time

similarities.df
