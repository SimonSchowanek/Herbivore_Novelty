######################################################################!
#                                                                    #
#      IAa6: SELECT NON-ANALOGUE CELLS                               #  
#                                                                    #
######################################################################!

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # METADATA
        # This script loads the earlier saved distance matrix and select which cells are non-analogue 
        #
        # Output: a table indicating which cells are non-analogue and a map visualising the results
        # time: needs 20 min to run.
        # Status: Ready (2025/02/03), but see notes
        # Note: This scrips used large matrices (2-10 Gb). This is heavier than R can work with. Therefore we use the Bigmemory package.
        #       This package requires loading the csv file distance matrix. This takes long the first time, but should be quicker afterwards.
        #
        #
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
library(geodata)


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

        #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
        # – Load Data ####
        #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

functional.community.matrix.relative = read.csv("./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_relative_abundances.csv")
functional.community.matrix.relative = functional.community.matrix.relative[,!(colnames(functional.community.matrix.relative) %in% c("X","CellNr"))] # remove two redundant variables

# Spatial
#–––––––––––––––––––––––––––––––––––––––––––––––––
raster.example = raster("./1_Input/PHYLACINE_V1_2_1/Ranges/Present_natural/Stegodon_orientalis.tif") # an example raster from PHYLACINE
biogeo.realms = readOGR("./1_Input/Biogeographical_Realms/udvardy/udvardy_py.shp")  # Shapefile of the biogeographical relams https://data-gis.unep-wcmc.org/portal/home/item.html?id=7f3a055f36104f36b6dd7834ebe2cf45
herbivomes = readOGR("./1_Input/Herbivome_Polygons.shp") # the herbivomes shapefile



## Herbivome Distances
df.herbivore.distances = read.csv("./3_Output/IAa5_Calculate_Distances_Herbivomes/IAa5a_Distance_Hempson_Herbivomes.csv")

# Attach (= load) the big distance matrix
#–––––––––––––––––––––––––––––––––––––––––––––––––
big.distance.matrix <- attach.big.matrix("./3_Output/IAa4_Calculate_Cell_Distances/IAa4_big_distance_matrix.desc")

# Check data
#–––––––––––––––––––––––––––––––––––––––––––––––––
big.distance.matrix[1:10,1:10]
dim(big.distance.matrix)
#range(big.distance.matrix[,]) #range = 0-2 (not a trivial calculation)

        #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
        # – Load Functions ####
        #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

source("./2_scripts/Functions/I_F1_Realm_Distances.R")

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

df.pn = df.master[which(df.master$Scenario == "Present-Natural"),]
df.current = df.master[df.master$Scenario == "Current",]



## – Directories  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
if(!dir.exists(paste0("./3_Output/IAa6_Map_Non_Analog_Communities"))){
  dir.create(paste0("./3_Output/IAa6_Map_Non_Analog_Communities"))}


#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 3. DIVIDE DISTANCE MATRIX ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# Divide cells between current and present-natural cells
Cells.Current = which(functional.community.matrix.relative[which(functional.community.matrix.relative$Herbivores.Present == 1),"Scenario"] == "Current")
Cells.PN = which(functional.community.matrix.relative[which(functional.community.matrix.relative$Herbivores.Present == 1),"Scenario"] == "Present-Natural")

length(Cells.PN) + length(Cells.Current)
        
                #~~~~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # 
                #   Squared chord distance values can range from 0.0 to 2.0, with 0.0 indicating identical proportions of species within the samples being compared. 
                #
                #   For the squared chord-distance metric, critical values have ranged from 0.05 in local applications to 0.40 for continental-scale paleoclimate reconstructions.
                #   values of 0.12 to 0.20 most commonly used in vegetational and paleoclimatic inference at intermediate scales.
                #   Several authors have used multiple critical values to indicate varying degrees of confidence in analog determination or to draw paleoenvironmental inferences at multiple ecological scales 
                #
                #    SOURCE:  Jackson, S. T., and J. W. Williams. 2004. MODERN ANALOGS IN QUATERNARY PALEOECOLOGY: Here Today, Gone Yesterday, Gone Tomorrow? Annual Review of Earth and Planetary Sciences 32:495–537.
                #
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
        # – Subset Data (distance matrix) ####
        #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# How similar are current cells to each other?
#–––––––––––––––––––––––––––––––––––––––––––––––––

Current.vs.Current <- sub.big.matrix(big.distance.matrix,
                                             firstRow = 1,
                                             lastRow = length(Cells.Current),
                                             firstCol = 1,
                                             lastCol = length(Cells.Current))
dim(Current.vs.Current)

# How similar are current cells to PN cells?
#–––––––––––––––––––––––––––––––––––––––––––––––––
Current.vs.PN <- sub.big.matrix(big.distance.matrix,
                                     firstRow = 1,
                                     lastRow = length(Cells.Current),
                                     firstCol = length(Cells.Current) + 1,
                                     lastCol = length(Cells.Current) + length(Cells.PN))
dim(Current.vs.PN)

# How similar are PN to other PN cells?
#–––––––––––––––––––––––––––––––––––––––––––––––––
PN.vs.PN <- sub.big.matrix(big.distance.matrix,
                                firstRow = length(Cells.Current) + 1,
                                lastRow = length(Cells.Current) + length(Cells.PN),
                                firstCol = length(Cells.Current) + 1,
                                lastCol = length(Cells.Current) + length(Cells.PN))

dim(PN.vs.PN)

# How similar are PN to current cells? 
#–––––––––––––––––––––––––––––––––––––––––––––––––
PN.vs.Current <- sub.big.matrix(big.distance.matrix,
                                firstRow = length(Cells.Current) + 1,
                                lastRow = length(Cells.Current) + length(Cells.PN),
                                firstCol = 1,
                                lastCol = length(Cells.Current))

dim(PN.vs.Current)

## Turn Diagonal into NA
#–––––––––––––––––––––––––––––––––––––––––––––––––
PN.vs.PN.with.NA = PN.vs.PN
d <- min(dim(PN.vs.PN.with.NA))
PN.vs.PN.with.NA[cbind(1:d, 1:d)] = NaN # Is the same as NA, but in C++

Current.vs.Current.with.NA= Current.vs.Current
d <- min(dim(Current.vs.Current.with.NA))
Current.vs.Current.with.NA[cbind(1:d, 1:d)] = NaN # Is the same as NA, but in C++

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 4. Map degree of functional change ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# Select cells that have herbivores in both current and pn
cells.in.current.and.pn = which(df.pn$Original.Cell %in% df.current$Original.Cell)
distance.between.current.and.pn = Current.vs.PN[,cells.in.current.and.pn]
dim(distance.between.current.and.pn)
distance.between.current.and.pn[1:10,1:10]

# Take the distance between the same cell in current and Present-natural
d <- min(dim(distance.between.current.and.pn))  # select smallest dimension
distance.between.current.and.pn = distance.between.current.and.pn[cbind(1:d, 1:d)]

# Attach the resulting distance to the df.current.df
df.current$Distance.from.pn = distance.between.current.and.pn

# Add Cells where herbivores have totally dissapeared
cells.in.pn.not.in.current = which(!(df.pn$Original.Cell %in% df.current$Original.Cell))
Areas.that.lost.all.herbivores = df.pn[cells.in.pn.not.in.current,]






#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 5. OCCURRENCE OF NON-ANALOGUE COMMUNITIES  ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

## Save table
df.info = data.frame(matrix(NA,1,8))
colnames(df.info) = c("threshold", "threshold.conservative", "non.analog.without.islands", "non.analog.with.islands", "perc.cells.changed", "perc.cells.lost.all.herbivores", "non.analog.with.islands.minus.0.1", "non.analog.with.islands.plus.0.1")

#define threshold
print(df.herbivore.distances)
df.info$threshold = min(df.herbivore.distances[which(df.herbivore.distances[,1] %in% c("Forest Duiker", "Arid Gazelle", "Savanna Combined")), "between"])
df.info$threshold.conservative = sort(df.herbivore.distances[which(df.herbivore.distances[,1] %in% c("Forest Duiker", "Arid Gazelle", "Savanna Combined")), "between"])[2]
df.info$threshold.liberal = 0.4


# Areas without novelty 
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > df.info$threshold),"threshold.value"] = 1


# How many cells are non-analogue (without islands) 
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > df.info$threshold),"threshold.value"] = 1
sum(df.current$threshold.value)/nrow(df.pn) #
round(sum(df.current$threshold.value)/length(df.current$threshold.value),2)
df.info[1,"non.analog.without.islands"] = round(sum(df.current$threshold.value)/length(df.current$threshold.value),2)


# How many cells are non-analogue (with islands) 
cell.that.lost.all.herbivores = length(which(!(df.pn$Original.Cell %in% df.current$Original.Cell)))
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > df.info$threshold),"threshold.value"] = 1
non.analogue.cells.perc = (sum(df.current$threshold.value) + cell.that.lost.all.herbivores)/nrow(df.pn) #
df.info[1,"non.analog.with.islands"] = round(non.analogue.cells.perc,2)


# How many cells are how seen some change?
cell.that.lost.all.herbivores = length(which(!(df.pn$Original.Cell %in% df.current$Original.Cell)))
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > 0),"threshold.value"] = 1
non.analogue.cells.perc = (sum(df.current$threshold.value) + cell.that.lost.all.herbivores)/nrow(df.pn) #
df.info[1,"perc.cells.changed"] = round(non.analogue.cells.perc,2)



# How many cells have a lost all their herbivores
cell.that.lost.all.herbivores = length(which(!(df.pn$Original.Cell %in% df.current$Original.Cell)))
cell.that.lost.all.herbivores.perc = cell.that.lost.all.herbivores/nrow(df.pn) #
df.info[1,"perc.cells.lost.all.herbivores"] = round(cell.that.lost.all.herbivores.perc,2) # cells that lost all herbivores
#round(sum(df.current$threshold.value)/length(df.current$threshold.value),2)


# How many cells are non-analogue (with islands) 
cell.that.lost.all.herbivores = length(which(!(df.pn$Original.Cell %in% df.current$Original.Cell)))
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > df.info$threshold - 0.1),"threshold.value"] = 1
non.analogue.cells.perc = (sum(df.current$threshold.value) + cell.that.lost.all.herbivores)/nrow(df.pn) #
df.info[1,"non.analog.with.islands.minus.0.1"] = round(non.analogue.cells.perc,2)



# How many cells are non-analogue (with islands) 
cell.that.lost.all.herbivores = length(which(!(df.pn$Original.Cell %in% df.current$Original.Cell)))
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > df.info$threshold + 0.1),"threshold.value"] = 1
non.analogue.cells.perc = (sum(df.current$threshold.value) + cell.that.lost.all.herbivores)/nrow(df.pn) #
df.info[1,"non.analog.with.islands.plus.0.1"] = round(non.analogue.cells.perc,2)







# How many cells are non-analogue (liberal and with islands) 
cell.that.lost.all.herbivores = length(which(!(df.pn$Original.Cell %in% df.current$Original.Cell)))
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > df.info$threshold.liberal), "threshold.value"] = 1
non.analogue.cells.perc = (sum(df.current$threshold.value) + cell.that.lost.all.herbivores)/nrow(df.pn) #
df.info[1,"non.analog.with.islands.liberal"] = round(non.analogue.cells.perc,2)



# How many cells are non-analogue (conservative and with islands) 
cell.that.lost.all.herbivores = length(which(!(df.pn$Original.Cell %in% df.current$Original.Cell)))
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > df.info$threshold.conservative),"threshold.value"] = 1
non.analogue.cells.perc = (sum(df.current$threshold.value) + cell.that.lost.all.herbivores)/nrow(df.pn) #
df.info[1,"non.analog.with.islands.conservative"] = round(non.analogue.cells.perc,2)


## Save info
write.csv(df.info, "./3_Output/IAa6_Map_Non_Analog_Communities/IAa6a_saved_parameters.csv")



#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 6. PLOT OCCURRENCE OF NON-ANALOGUE COMMUNITIES  ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

## CONTINUOUS
#–––––––––––––––––––––––––––––––––––––––––––––––––

colours.1 = brewer.pal(4, "PuBu")
colours.2 = brewer.pal(5, "OrRd")
colours.final = c(colours.1,colours.2)
pie(rep(1,length(colours.final)), col=colours.final)



## DISCRETE 2025 NORMAL VERSION 
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Turn the continuous dissimilarity values into bins
threshold = df.info$threshold 


## Breaks
(breaks = sort(c(seq(0,2, by = 0.25), threshold)))
length(breaks)


## Turn data into bins
df.current$Distance.Discrete = cut(df.current$Distance.from.pn,
                                   breaks = breaks, 
                                   right = F,
                                   include.lowest = T)
str(df.current)


## Create a colour scheme 
nbins = length(unique(levels(df.current$Distance.Discrete)))
n.bins.below = length(which(breaks<=threshold))-1 # nr of bins below the treshold
n.bins.above = length(which(breaks>=threshold))-1 

colours.1 = brewer.pal(nbins, "PuBu") # blues
colours.1 = colours.1[3:(2 + n.bins.below)]
colours.2 = brewer.pal(nbins, "OrRd") # reds
colours.2 = colours.2[(n.bins.below+1):nbins]

colours.final = c(colours.1,colours.2)
pie(rep(1,length(colours.final)), col=colours.final)


## Plot
latitude.line = -2.8*10^5 # sweet spot for pdf
latitude.line = -7.5*10^5# sweet spot for TIFF
    
(p1 = ggplot(data=df.current, aes(y=Latitude, x=Longitude, fill = Distance.Discrete), color = Distance.Discrete) +
    
    geom_tile(size=0.5) +
    
    scale_fill_manual(values = colours.final,
                      drop = F) +
    scale_colour_manual(values = colours.final) +
    
    geom_polygon(data = countries.poly, 
                 aes(x = long, 
                     y = lat, 
                     group = group),
                 fill = NA, 
                 colour = 'gray10', 
                 size = 0.3) +
  
    coord_fixed() +
    theme_void()  +
    
    #    ggtitle("Current VS present-natural") + # add title
    labs(x="", y="", fill = "Dissimilarity", title = paste0("Distribution of novel herbivore assemblages \n (treshold = ", df.info$threshold, ")")) +
    
    theme(
          plot.title = element_text("sans", "plain", "Black", 24, hjust = 0.5),
          legend.title = element_blank(),
          legend.key = element_rect(color = NA),
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.text = element_text(size=20),
          legend.position = c(0.1, 0.35),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.background = element_blank()) +
    
    geom_tile(data=Areas.that.lost.all.herbivores, aes(y=Latitude, x=Longitude), inherit.aes = F, colour ="black", fill = "black"))
    

## Save PDF
filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_discrete_threshold_", threshold)
## Save PDF
pdf(width=15,height=8, file = paste0(filename, ".pdf"))
p1
dev.off()

## Save PNG
png(width=1400,height=700, file = paste0(filename, ".png"))
p1
dev.off()







## DISCRETE 2025 CONSERVATIVE VERSION 
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Turn the continuous dissimilarity values into bins
threshold = df.info$threshold.conservative

breaks = sort(c(seq(0,2, by = 0.25), threshold))
length(breaks)

df.current$Distance.Discrete = cut(df.current$Distance.from.pn,
                                   breaks = breaks, 
                                   right = F,
                                   include.lowest = T)
str(df.current)



## Create a colour scheme 
nbins = length(unique(levels(df.current$Distance.Discrete)))
n.bins.below = length(which(breaks<=threshold))-1 # nr of bins below the treshold
n.bins.above = length(which(breaks>=threshold))-1 

colours.1 = brewer.pal(nbins, "PuBu") # blues
colours.1 = colours.1[3:(2 + n.bins.below)]
colours.2 = brewer.pal(nbins, "OrRd") # reds
colours.2 = colours.2[(n.bins.below+1):nbins]

colours.final = c(colours.1,colours.2)
pie(rep(1,length(colours.final)), col=colours.final)



## Plot
latitude.line = -2.8*10^5 # sweet spot for pdf
latitude.line = -7.5*10^5# sweet spot for TIFF

(p2 = ggplot(data=df.current, aes(y=Latitude, x=Longitude, fill = Distance.Discrete), color = Distance.Discrete) +
    geom_tile(size=0.5) +
    scale_fill_manual(values = colours.final,
                      drop = F) +
    scale_colour_manual(values = colours.final) +
    geom_polygon(data = countries.poly, 
                 aes(x = long, 
                     y = lat, 
                     group = group),
                 fill = NA, 
                 colour = 'gray10', 
                 size = 0.3) +
    
    
    coord_fixed() +
    
    theme_void()  +
    
    #    ggtitle("Current VS present-natural") + # add title
    labs(x="", y="", fill = "Dissimilarity", title = paste0("Distribution of novel herbivore assemblages \n (treshold = ", df.info$threshold.conservative, ")")) +
    
    theme(plot.title = element_text("sans", "plain", "Black", 24, hjust = 0.5),
          legend.title = element_blank(),
          legend.key = element_rect(color = NA),
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.text = element_text(size=20),
          legend.position = c(0.1, 0.35),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.background = element_blank()) +
    geom_tile(data=Areas.that.lost.all.herbivores, aes(y=Latitude, x=Longitude), inherit.aes = F, colour ="black", fill = "black"))
    # geom_segment(aes(x = (min(Longitude)-5*10^5), xend = (min(Longitude) + 9*10^6), y = latitude.line, yend = latitude.line), colour = "black", size = 1) +
    # annotate("text", x = (-1.15*10^7), y = latitude.line + 1*10^6, label = "Analog", size = 10) +
    # annotate("text", x = (-1.15*10^7), y = latitude.line - 1*10^6, label = "Non-Analog", size = 10))

## Save PDF
filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_discrete_threshold_", threshold)
## Save PDF
pdf(width=15,height=8, file = paste0(filename, ".pdf"))
p2
dev.off()

## Save PNG
png(width=1400,height=700, file = paste0(filename, ".png"))
p2
dev.off()










## Save PNG
filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_discrete_threshold_combined")
png(width=1000,height=1000, file = paste0(filename, ".png"))
ggarrange(p1,
          p2,
          nrow = 2, 
          labels = "AUTO")
dev.off()

















## DISCRETE 2025 LIBERAL VERSION 
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Turn the continuous dissimilarity values into bins
df.info$threshold.liberal= 0.4
threshold = df.info$threshold.liberal

breaks = sort(c(seq(0,2, by = 0.25), threshold))
length(breaks)

df.current$Distance.Discrete = cut(df.current$Distance.from.pn,
                                   breaks = breaks, 
                                   right = F,
                                   include.lowest = T)
str(df.current)



## Create a colour scheme 
nbins = length(unique(levels(df.current$Distance.Discrete)))
n.bins.below = length(which(breaks<=threshold))-1 # nr of bins below the treshold
n.bins.above = length(which(breaks>=threshold))-1 

colours.1 = brewer.pal(nbins, "PuBu") # blues
colours.1 = colours.1[3:(2 + n.bins.below)]
colours.2 = brewer.pal(nbins, "OrRd") # reds
colours.2 = colours.2[(n.bins.below+1):nbins]

colours.final = c(colours.1,colours.2)
pie(rep(1,length(colours.final)), col=colours.final)



## Plot
latitude.line = -2.8*10^5 # sweet spot for pdf
latitude.line = -7.5*10^5# sweet spot for TIFF

(p3 = ggplot(data=df.current, aes(y=Latitude, x=Longitude, fill = Distance.Discrete), color = Distance.Discrete) +
    geom_tile(size=0.5) +
    scale_fill_manual(values = colours.final,
                      drop = F) +
    scale_colour_manual(values = colours.final) +
    geom_polygon(data = countries.poly, 
                 aes(x = long, 
                     y = lat, 
                     group = group),
                 fill = NA, 
                 colour = 'gray10', 
                 size = 0.3) +
    
    
    coord_fixed() +
    
    theme_void()  +
    
    #ggtitle("Distribution of novel herbivore assemblages") + # add title
    
    labs(x="", y="", fill = "Dissimilarity", title = paste0("Distribution of novel herbivore assemblages \n (treshold = ", df.info$threshold.liberal, ")")) +
    
    theme(plot.title = element_text("sans", "plain", "Black", 24, hjust = 0.5),
          legend.title = element_blank(),
          legend.key = element_rect(color = NA),
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.text = element_text(size=20),
          legend.position = c(0.1, 0.35),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.background = element_blank()) +
    geom_tile(data=Areas.that.lost.all.herbivores, aes(y=Latitude, x=Longitude), inherit.aes = F, colour ="black", fill = "black"))
# geom_segment(aes(x = (min(Longitude)-5*10^5), xend = (min(Longitude) + 9*10^6), y = latitude.line, yend = latitude.line), colour = "black", size = 1) +
# annotate("text", x = (-1.15*10^7), y = latitude.line + 1*10^6, label = "Analog", size = 10) +
# annotate("text", x = (-1.15*10^7), y = latitude.line - 1*10^6, label = "Non-Analog", size = 10))

## Save PDF
filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_discrete_threshold_", threshold)
## Save PDF
pdf(width=15,height=8, file = paste0(filename, ".pdf"))
p3
dev.off()

## Save PNG
png(width=1400,height=700, file = paste0(filename, ".png"))
p3
dev.off()










## Save PNG
filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_discrete_threshold_combined")
png(width=1000,height=1000, file = paste0(filename, ".png"))
ggarrange(p3,
          p2,
          nrow = 2, 
          labels = "AUTO")
dev.off()














#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 10. CREATE MAP WITH UNCERTAIN ECOSYSTEMS ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# – Load Data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––


    #~~~~~ NOTE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #       
    #   2025/02/03: This section used the Raster package, which no longer seems to work. I have adapted the script so that it runs again.
    #               The end result should be identical but I am leaving this note here for future reference
    # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create directory to save data
if(!dir.exists(paste0("3_Output/IAa6_Map_Non_Analog_Communities/climate_data"))){ # create directory
  dir.create(paste0("3_Output/IAa6_Map_Non_Analog_Communities/climate_data"))}


# Download WorldClim bioclimatic data at 10-minute resolution (~18 km)
r <- worldclim_global(var = "bio", res = 10, path = "3_Output/IAa6_Map_Non_Analog_Communities/climate_data")


#Bio 1 and Bio12 are mean anual temperature and anual precipitation:(https://wec.wur.nl/r/spatial/raster-data.html)
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec") # give names.


## – Set parameters ####
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Visual Parameters
theme_set(theme_bw()) # set ggplot layout

## – Correct data ####
#–––––––––––––––––––––––––––––––––––––––––––––––––

# Ecosystems uncertain were mapped for all pixels with MAP > 7.143 MAT + 286 and MAP < –1.469 MAT2 + 81.665 MAT + 475
Con1 = r$Prec > 7.143 * r$Temp + 286
Con2 = r$Prec < -1.469 * r$Temp^2 + 81.665 * r$Temp + 475

plot(Con1)
plot(Con2)
r.combined = Con1 + Con2
r.combined = r.combined >= 2 # select areas that meet both conditions
plot(r.combined) # Uncertain Areas


library(terra)
# Reproject r.combined to match raster.example
r.uncertain.repr <- project(r.combined, raster.example)

# Check the output
print(r.uncertain.repr)
plot(r.uncertain.repr[[1]])  # Plot first layer

r.uncertain.repr = r.uncertain.repr > 0 # include uncertain areas
plot(r.uncertain.repr)



# Add variable to indicate which cells are uncertain
v.uncertain.repr = as.numeric(values(r.uncertain.repr))
df.current$Uncertain = v.uncertain.repr[df.current$Original.Cell]

# make changes so that area that lost all herbivores and had uncertain ecosystems are also included.
df.uncertain.repr = data.frame(1:length(v.uncertain.repr), v.uncertain.repr)
matches = which(Areas.that.lost.all.herbivores$Original.Cell %in% df.uncertain.repr[which(df.uncertain.repr$v.uncertain.repr == T), "X1.length.v.uncertain.repr."]) # select cells that have lost all herbivores and that are uncertain ecosystems
Areas.that.lost.all.herbivores.and.were.uncertain = Areas.that.lost.all.herbivores[matches,]

## – Plot data ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
df.current[which((df.current$Distance.from.pn >= df.info$threshold) & (df.current$Uncertain == T)),"Status"] = "ABS possible, large herbivore change"
df.current[which((df.current$Distance.from.pn < df.info$threshold) & (df.current$Uncertain == T)),"Status"] = "ABS possible, small herbivore change" 
df.current[which((df.current$Distance.from.pn >= df.info$threshold) & (df.current$Uncertain == F)),"Status"] = "ABS impossible, large herbivore change" 



(df.info$novel.herbivory.abs.possible = length(which(df.current$Status== "Alternative biome States possible, large herbivore change")) / nrow(df.current))
(df.info$no.novel.herbivory.abs.possible = length(which(df.current$Status== "Alternative biome States possible, small herbivore change")) / nrow(df.current))
(df.info$novel.herbivory.abs.impossible = length(which(df.current$Status== "Alternative biome States impossible, large herbivore change")) / nrow(df.current))


## Select sites with 
df.current.subset = df.current[!is.na(df.current$Status),]

(p.veg.change = ggplot(data=df.current.subset, aes(y=Latitude, x=Longitude, fill = Status)) +
    geom_tile(size=0.5) +
    scale_fill_manual(values=c("#FFB000", "#648FFF", "#DC267F")) +
    geom_polygon(data = countries.poly, 
                 aes(x = long, 
                     y = lat, 
                     group = group),
                 fill = NA, 
                 colour = 'gray10', 
                 size = 0.3) +
    coord_fixed() +
    theme_void()  +
    ggtitle("Predicted vegetation response to novel herbivore impacts") + # add title
    labs(x="", y="", fill = "Dissimilarity") +
    theme(plot.title = element_text("sans", "plain", "Black", 18, hjust = 0.5),
          legend.position = "bottom",
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.background = element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(size=20)) +
    geom_tile(data=Areas.that.lost.all.herbivores.and.were.uncertain, aes(y=Latitude, x=Longitude), inherit.aes = F, colour ="#648FFF", fill = "#648FFF"))
#+                guides(fill=guide_legend(nrow=3, byrow=TRUE))

## Save PDF
filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_uncertain_with_change")
## Save PDF
pdf(width=15,height=8, file = paste0(filename, ".pdf"))
p.veg.change
dev.off()

## Save PNG
png(width=1400,height=700, file = paste0(filename, ".png"))
p.veg.change
dev.off()


## Combine with Non-analog plot

## Save PNG
filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_uncertain_with_change_combined")
png(width=1400,height=1400, file = paste0(filename, ".png"))
ggarrange(p1, p.veg.change, nrow = 2, labels = "AUTO")
dev.off()

filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_uncertain_with_change_combined")
pdf(width=14,height=14, file = paste0(filename, ".pdf"))
ggarrange(p1, p.veg.change, nrow = 2, labels = "AUTO")
dev.off()






## – End Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
end_time <- Sys.time()
end_time - start_time

