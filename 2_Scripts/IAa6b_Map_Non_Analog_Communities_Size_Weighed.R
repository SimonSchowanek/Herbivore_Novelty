######################################################################!
#                                                                    #
#      IAa6: SELECT NON-ANALOGUE CELLS     (SIZE WEIGHED)            #  
#                                                                    #
######################################################################!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# METADATA
# This script loads the earlier saved distance matrix and select which cells are non-analogue 
#
# Output: a table indicating which cells are non-analogue and a map visualising the results
# time: needs 2 min to run.
# Status: Ready, Dasypus bellis has been removed (23/01/2020)
# Note: This scrips used large matrices (2-10 Gb). This is heavier than R can work with. Therefore we use the Bigmemory package.
#       This package requires loading the csv file distance matrix. This takes long the first time, but should be quicker afterwards.
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
library(ggpubr) # multiple panels in one plot


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

functional.community.matrix.relative = read.csv("./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_relative_abundances_weighed.csv")
functional.community.matrix.relative = functional.community.matrix.relative[,!(colnames(functional.community.matrix.relative) %in% c("X","CellNr"))] # remove two redundant variables

## The species community matrix
community.matrix = read.csv("./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_community_matrix_combined_species.csv")
colnames(community.matrix)

# Spatial
#–––––––––––––––––––––––––––––––––––––––––––––––––
raster.example = raster("/Users/au572919/Datasets/PHYLACINE_V1_2_1/Ranges/Current/Abditomys_latidens.tif")
biogeo.realms = readOGR("/Users/au572919/Datasets/Biogeographical_Realms/udvardy/udvardy_py.shp")
herbivomes = readOGR("./1_files/Herbivome_Polygons.shp")

## Herbivome Distances
df.herbivore.distances = read.csv("./3_Output/IAa5_Calculate_Distances_Herbivomes/IAa5a_Distance_Hempson_Herbivomes_weighed.csv")

# Attach (= load) the big distance matrix
#–––––––––––––––––––––––––––––––––––––––––––––––––
big.distance.matrix <- attach.big.matrix("./3_Output/IAa4_Calculate_Cell_Distances/IAa4_big_distance_matrix_weighed.desc")

# Check data
#–––––––––––––––––––––––––––––––––––––––––––––––––
big.distance.matrix[1:10,1:10]
dim(big.distance.matrix)
#range(big.distance.matrix[,]) #range = 0-2 (not a trivial calculation)

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


## – Set parameters ####
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Visual Parameters
theme_set(theme_bw()) # set ggplot layout

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
df.info = data.frame(matrix(NA,1,7))
colnames(df.info) = c("threshold", "non.analog.without.islands", "non.analog.with.islands", "perc.cells.changed", "perc.cells.lost.all.herbivores", "non.analog.with.islands.minus.0.1", "non.analog.with.islands.plus.0.1")

#define threshold
print(df.herbivore.distances)
df.info$threshold = min(df.herbivore.distances[which(df.herbivore.distances[,1] %in% c("Forest Duiker", "Arid Gazelle", "Savanna Combined")), "between"])

# Areas without change 
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


#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 6. PLOT OCCURRENCE OF NON-ANALOGUE COMMUNITIES  ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

## CONTINUOUS
#–––––––––––––––––––––––––––––––––––––––––––––––––

colours.1 = brewer.pal(4, "PuBu")
colours.2 = brewer.pal(5, "OrRd")
colours.final = c(colours.1,colours.2)
pie(rep(1,length(colours.final)), col=colours.final)

#threshold = 0.55

# pdf(width=15,height=8, file = "3_Output/IAa6_Map_Non_analogue_communities/IAa6b_non_analogues_communities_world_threshold_70.pdf")
# 
# ggplot(data=df.current, aes(y=Latitude, x=Longitude, fill = Distance.from.pn), color = Distance.from.pn) +
#     geom_tile(size=0.5) +
#     scale_fill_gradientn(limits = c(0, 2), colours=colours.final,
#                         breaks = c(0,threshold,1,2),
#                          values= rescale(c(0,(threshold)-1*10^(-1.4),threshold, 2))) +
#     # scale_colour_gradientn(limits = c(0, 2), colours=colours.final,
#     #                       breaks = c(0,threshold,1,2),v
#     #                       values= rescale(c(0,threshold-1*10^(-10),threshold, 2))) +
#     geom_polygon(data = countries.poly, 
#                  aes(x = long, 
#                      y = lat, 
#                      group = group),
#                  fill = NA, 
#                  colour = 'gray10', 
#                  size = 0.3) +
#     coord_fixed() +
#     theme_void()  +
#     ggtitle("Current VS present-natural") + # add title
#     labs(x="", y="", fill = "Dissimilarity") +
#     theme(plot.title = element_text("sans", "plain", "Black", 18, hjust = 0.5),
#           legend.title = element_blank(),
#           legend.position="bottom",
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           legend.background = element_rect(color = "white", fill = "white", size = 0.5, linetype = "solid")) +
#     guides(fill = guide_colourbar(barwidth = 20, barheight = 1)) +
#     geom_tile(data=Areas.that.lost.all.herbivores, aes(y=Latitude, x=Longitude), inherit.aes = F, colour ="black", fill = "black") # cells where herbivores disappeared
# dev.off()


## DISCRETE 2021 VERSION (September)
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Turn the continous dissimilarity values into bins

# breaks = c(0,
#            threshold*1/3,
#            threshold*2/3,
#            threshold,
#            threshold + (2-threshold)*1/6,
#            threshold + (2-threshold)*2/6,
#            threshold + (2-threshold)*3/6,
#            threshold + (2-threshold)*4/6,
#            threshold + (2-threshold)*5/6,
#            2)

#threshold = 0.55

breaks = sort(c(seq(0,2, by = 0.25), df.info$threshold))
length(breaks)

df.current$Distance.Discrete = cut(df.current$Distance.from.pn,
                                   breaks = breaks, 
                                   right = F,
                                   include.lowest = T)
str(df.current)

## Create a colour scheme 
colours.1 = brewer.pal(9, "PuBu")
colours.1 = colours.1[3:5]
colours.2 = brewer.pal(9, "OrRd")
colours.2 = colours.2[4:9]

colours.final = c(colours.1,colours.2)
pie(rep(1,length(colours.final)), col=colours.final)


## Plot
latitude.line = -2.8*10^5 # sweet spot for pdf
latitude.line = -7.5*10^5# sweet spot for TIFF

(p.functional.change = ggplot(data=df.current, aes(y=Latitude, x=Longitude, fill = Distance.Discrete), color = Distance.Discrete) +
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
        labs(x="", y="", fill = "Dissimilarity") +
        theme(plot.title = element_text("sans", "plain", "Black", 18, hjust = 0.5),
              legend.title = element_blank(),
              legend.key = element_rect(color = NA),
              legend.key.height = unit(1.5, 'cm'), #change legend key height
              legend.key.width = unit(1.5, 'cm'), #change legend key width
              legend.text = element_text(size=16),
              legend.position = c(0.1, 0.3),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.background = element_blank()) +
        geom_tile(data=Areas.that.lost.all.herbivores, aes(y=Latitude, x=Longitude), inherit.aes = F, colour ="black", fill = "black") +
        geom_segment(aes(x = (min(Longitude)-5*10^5), xend = (min(Longitude) + 9*10^6), y = latitude.line, yend = latitude.line), colour = "black", size = 1) +
        annotate("text", x = (-1.15*10^7), y = latitude.line + 1*10^6, label = "Analog", size = 10) +
        annotate("text", x = (-1.15*10^7), y = latitude.line - 1*10^6, label = "Non-Analog", size = 10))

## Save PDF
filename = paste0("3_Output/IAa6_Map_Non_Analog_Communities/IAa6_non_analogues_communities_world_discrete_threshold_weighed")
## Save PDF
pdf(width=15,height=8, file = paste0(filename, ".pdf"))
p.functional.change
dev.off()

## Save PNG
png(width=1400,height=700, file = paste0(filename, ".png"))
p.functional.change
dev.off()

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 7. GLOBALLY NON-ANALOGUE  ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

PN.vs.PN.subset.normal.matrix = as.matrix(PN.vs.PN)
Current.vs.PN.subset.normal.matrix = as.matrix(Current.vs.PN)

## PN VS PN
system.time(Quantile.PN.vs.PN.95.perc <- quantile(PN.vs.PN.subset.normal.matrix,0.95, na.rm = T))
Quantile.PN.vs.PN.95.perc

system.time(Quantile.PN.vs.PN.90.perc <- quantile(PN.vs.PN.subset.normal.matrix,0.90, na.rm = T))
Quantile.PN.vs.PN.90.perc

system.time(Quantile.PN.vs.PN.69.perc <- quantile(PN.vs.PN.subset.normal.matrix,0.69, na.rm = T))
Quantile.PN.vs.PN.69.perc # 69% of all distances is smaller than 0.97%

which(Current.vs.PN > Quantile.PN.vs.PN.69.perc)

## Current VS PN
system.time(Quantile.Current.vs.PN.95.perc <- quantile(Current.vs.PN.subset.normal.matrix,0.95, na.rm = T))
Quantile.Current.vs.PN.95.perc

system.time(Quantile.Current.vs.PN.90.perc <- quantile(Current.vs.PN.subset.normal.matrix,0.90, na.rm = T))
Quantile.Current.vs.PN.90.perc

system.time(Quantile.Current.vs.PN.69.perc <- quantile(Current.vs.PN.subset.normal.matrix,0.69, na.rm = T))
Quantile.Current.vs.PN.69.perc

#~~~~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Within all PN cells, 95% of the distances are smaller than 1.563564 ( = Quantile.PN.vs.PN.95.perc). 
# If Current.VS.PN are even more "distant" it means they are effectively "Globally Novel".        
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# How many cells are non-analogue (without islands) 
df.current[,"threshold.value"] = df.info$threshold
df.current[which(df.current$Distance.from.pn > Quantile.PN.vs.PN.95.perc),"threshold.value"] = 1
sum(df.current$threshold.value)/nrow(df.pn) #
round(sum(df.current$threshold.value)/length(df.current$threshold.value),2)

# How many cells are non-analogue (with islands) 
cell.that.lost.all.herbivores = length(which(!(df.pn$Original.Cell %in% df.current$Original.Cell)))
df.current[,"threshold.value"] = 0
df.current[which(df.current$Distance.from.pn > Quantile.PN.vs.PN.90.perc),"threshold.value"] = 1
non.analogue.cells.perc = (sum(df.current$threshold.value) + cell.that.lost.all.herbivores)/nrow(df.pn) #
round(non.analogue.cells.perc,2)    


#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 8. REGIONAL NON-ANALOGUE  ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####


## Present-Natural Distances for all biogeographical realms
#–––––––––––––––––––––––––––––––––––––––––––––––––

df.quantiles.realms = data.frame(matrix(NA,6,3))
colnames(df.quantiles.realms) = c("realm", "threshold.69%", "threshold.95%")
df.quantiles.realms$realm = c("Africotropical", "Indomalayan", "Palaearctic", "Nearctic", "Neotropical", "Australasian")

#~~~~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#   This section shows how heterogeneous PN assemblages would be (within the same biogeograhical realm)
#   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

matches = which(df.pn$Biogeo.realm == "Africotropical")
dist.PN.Africotropical = PN.vs.PN.with.NA[matches,matches]
hist(dist.PN.Africotropical)
df.quantiles.realms[which((df.quantiles.realms$realm == "Africotropical")), "threshold.69%"] = quantile(dist.PN.Africotropical,0.69, na.rm = T)
df.quantiles.realms[which((df.quantiles.realms$realm == "Africotropical")), "threshold.95%"] = quantile(dist.PN.Africotropical,0.95, na.rm = T)

matches = which(df.pn$Biogeo.realm == "Indomalayan")
dist.PN.Indomalayan = PN.vs.PN.with.NA[matches,matches]
hist(dist.PN.Indomalayan)
df.quantiles.realms[which((df.quantiles.realms$realm == "Indomalayan")), "threshold.69%"] = quantile(dist.PN.Indomalayan,0.69, na.rm = T)
df.quantiles.realms[which((df.quantiles.realms$realm == "Indomalayan")), "threshold.95%"] = quantile(dist.PN.Indomalayan,0.95, na.rm = T)

matches = which(df.pn$Biogeo.realm == "Palaearctic")
dist.PN.Palaearctic = PN.vs.PN.with.NA[matches,matches]
hist(dist.PN.Palaearctic)
df.quantiles.realms[which((df.quantiles.realms$realm == "Palaearctic")), "threshold.69%"] = quantile(dist.PN.Palaearctic,0.69, na.rm = T)
df.quantiles.realms[which((df.quantiles.realms$realm == "Palaearctic")), "threshold.95%"] = quantile(dist.PN.Palaearctic,0.95, na.rm = T)

matches = which(df.pn$Biogeo.realm == "Nearctic")
dist.PN.Nearctic = PN.vs.PN.with.NA[matches,matches]
hist(dist.PN.Nearctic)
df.quantiles.realms[which((df.quantiles.realms$realm == "Nearctic")), "threshold.69%"] = quantile(dist.PN.Nearctic,0.69, na.rm = T)
df.quantiles.realms[which((df.quantiles.realms$realm == "Nearctic")), "threshold.95%"] = quantile(dist.PN.Nearctic,0.95, na.rm = T)


matches = which(df.pn$Biogeo.realm == "Neotropical")
dist.PN.Neotropical = PN.vs.PN.with.NA[matches,matches]
hist(dist.PN.Neotropical)
df.quantiles.realms[which((df.quantiles.realms$realm == "Neotropical")), "threshold.69%"] = quantile(dist.PN.Neotropical,0.69, na.rm = T)
df.quantiles.realms[which((df.quantiles.realms$realm == "Neotropical")), "threshold.95%"] = quantile(dist.PN.Neotropical,0.95, na.rm = T)


matches = which(df.pn$Biogeo.realm == "Australasian")
dist.PN.Australasian = PN.vs.PN.with.NA[matches,matches]
hist(dist.PN.Australasian)
df.quantiles.realms[which((df.quantiles.realms$realm == "Australasian")), "threshold.69%"] = quantile(dist.PN.Australasian,0.69, na.rm = T)
df.quantiles.realms[which((df.quantiles.realms$realm == "Australasian")), "threshold.95%"] = quantile(dist.PN.Australasian,0.95, na.rm = T)

## Current Distances vs Present-Natural -  for all biogeographical realms
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Create New Variables
df.current[,"regionally.novel"] = 0

## Africotropical
realm = "Africotropical" # Select Realm
df.temp = df.current[which(df.current$Biogeo.realm == realm),] # Make new Df
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.69%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 1
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.95%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 2
df.current[which(df.current$Biogeo.realm == realm),"regionally.novel"] = df.temp$regionally.novel ## Overwrite the original df

## Indomalayan
realm = "Indomalayan" # Select Realm
df.temp = df.current[which(df.current$Biogeo.realm == realm),] # Make new Df
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.69%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 1
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.95%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 2
df.current[which(df.current$Biogeo.realm == realm),"regionally.novel"] = df.temp$regionally.novel ## Overwrite the original df

## Palaearctic
realm = "Palaearctic" # Select Realm
df.temp = df.current[which(df.current$Biogeo.realm == realm),] # Make new Df
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.69%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 1
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.95%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 2
df.current[which(df.current$Biogeo.realm == realm),"regionally.novel"] = df.temp$regionally.novel ## Overwrite the original df

## Nearctic
realm = "Nearctic" # Select Realm
df.temp = df.current[which(df.current$Biogeo.realm == realm),] # Make new Df
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.69%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 1
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.95%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 2
df.current[which(df.current$Biogeo.realm == realm),"regionally.novel"] = df.temp$regionally.novel ## Overwrite the original df

## Neotropical
realm = "Neotropical" # Select Realm
df.temp = df.current[which(df.current$Biogeo.realm == realm),] # Make new Df
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.69%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 1
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.95%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 2
df.current[which(df.current$Biogeo.realm == realm),"regionally.novel"] = df.temp$regionally.novel ## Overwrite the original df

## Australasia
realm = "Australasian" # Select Realm
df.temp = df.current[which(df.current$Biogeo.realm == realm),] # Make new Df
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.69%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 1
matches = which(df.temp[, "Distance.from.pn"] >= df.quantiles.realms[which(df.quantiles.realms$realm == realm), "threshold.95%"]) ## Identify Regionally Novel Assemblages
df.temp[matches, "regionally.novel"] = 2
df.current[which(df.current$Biogeo.realm == realm),"regionally.novel"] = df.temp$regionally.novel ## Overwrite the original df


## Current Distances vs Present-Natural -  for all biogeographical realms
#–––––––––––––––––––––––––––––––––––––––––––––––––

colours.final = c("#999999", "#E69F00", "firebrick")

(p2 = ggplot(data=df.current, aes(y=Latitude, x=Longitude, fill = as.factor(regionally.novel)), color =  as.factor(regionally.novel)) +
        geom_tile(size=0.5) +
        scale_fill_manual(values = colours.final,
                          drop = F,
                          labels = c("analog", "non-analog (69%)", "non-analog (95%)")) +
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
        labs(x="", y="", fill = "Dissimilarity") +
        theme(plot.title = element_text("sans", "plain", "Black", 18, hjust = 0.5),
              legend.title = element_blank(),
              legend.key = element_rect(color = NA),
              legend.key.height = unit(1.5, 'cm'), #change legend key height
              legend.key.width = unit(1.5, 'cm'), #change legend key width
              legend.text = element_text(size=15),
              legend.position = c(0.1, 0.3),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.background = element_blank()) +
        geom_tile(data=Areas.that.lost.all.herbivores, aes(y=Latitude, x=Longitude), inherit.aes = F, colour ="black", fill = "black"))


## Save TIFF
png(width=1400,height=700, file = "./3_Output/IAa6_Map_Non_Analog_Communities/IAa6b_non_analogues_communities_regionally_novel_weighed.png")
p2
dev.off()



## Plot differences
hist(Current.vs.PN.subset.normal.matrix, col = rgb(0,0,1,0.5), xlab = "Squared Chord Distance", main = "Distances 'PN vs PN' and 'Current vs PN'") # blue
hist(PN.vs.PN.subset.normal.matrix, col = rgb(1,0,0,0.5), add = T) # red
abline(v = 0.55, col = "black", add = T, lwd = 4, lty = "dashed")

## Save info
df.info$regionally.non.analog.69.perc = length(which(df.current$regionally.novel >= 1)) / nrow(df.current)
df.info$regionally.non.analog.95.perc = length(which(df.current$regionally.novel >= 2)) / nrow(df.current)


#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 9. COMPARE SPECIES RICHNESS LOSSES WITH FUNCTIONAL CHANGES  ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####


## Calculate species richness
columns.with.species = which(grepl("_", colnames(community.matrix))) # select all colnames with an underscore (only species names have these)
community.matrix$species.sum = apply(community.matrix[,columns.with.species],1,sum, na.rm = T)

# Remove species columns
community.matrix = community.matrix[,-columns.with.species]

# Divide matrix into current and present-natural
community.matrix.current = community.matrix[which(community.matrix$Scenario == "Current"),]
community.matrix.PN = community.matrix[which(community.matrix$Scenario == "Present-Natural"),]

## calculate species losses in current cells
community.matrix.current$pn.species.lost = community.matrix.PN$species.sum - community.matrix.current$species.sum
community.matrix.current$pn.species.lost.relative = 100 * (community.matrix.PN$species.sum - community.matrix.current$species.sum) / community.matrix.PN$species.sum

## Select cells where herbivores are present
community.matrix.current = community.matrix.current[community.matrix.current$Herbivores.Present == 1,]

## Functional Change
community.matrix.current$Distance.from.pn = df.current$Distance.from.pn

# community.matrix.current.2 = community.matrix.current[-which(community.matrix.current$pn.species.lost == 0),]

samples.change.tax.graph = sample(1:nrow(community.matrix.current), size = 2500)
community.matrix.current.2 = community.matrix.current[samples.change.tax.graph,]



## Plot Functional Change vs SR Loss
(p5 = ggplot(data = community.matrix.current.2, aes(x=pn.species.lost, y=Distance.from.pn)) +
        geom_point(aes(alpha = 0.005)) +
        geom_segment(aes(x = 0, y = 0, xend = max(pn.species.lost), yend =  max(Distance.from.pn), colour = "red")) +
        scale_x_continuous("Taxonomic Change (Species Richness Loss)") +
        scale_y_continuous("Functional Change (Squared Chord Distance)", limits = c(0,2)) +
        theme(legend.position='none') +
        annotate("label", x = 26, y = 0.1, label = "species richness loss overestimates functional change", size = 5, colour = "red") +
        annotate("label", x = 14, y = 1.7, label = "species richness loss underestimates functional change", size = 5, colour = "red")
)

(p6 = ggplot(data = community.matrix.current.2, aes(x=pn.species.lost.relative, y=Distance.from.pn)) +
                geom_point(aes(alpha = 0.005)) +
                geom_segment(aes(x = 0, y = 0, xend = 100, yend =  2, colour = "red")) +
                scale_x_continuous("Taxonomic Change (% of assemblage lost)") +
                scale_y_continuous("Functional Change (Squared Chord Distance)") +
                theme(legend.position='none') +
                annotate("label", x = 50, y = 0.1, label = "species richness loss overestimates functional change", size = 5, colour = "red") +
                annotate("label", x = 50, y = 1.7, label = "species richness loss underestimates functional change", size = 5, colour = "red")
)



## Save TIFF

png(width=600,height=600, file = "./3_Output/IAa6_Map_Non_Analog_Communities/IAa6b_species_loss_vs_functional_loss.png")
p6
dev.off()


## Save Info

df.info$nr.samples.function.taxonomy.graph = length(samples.change.tax.graph)
write.csv(df.info, "./3_Output/IAa6_Map_Non_Analog_Communities/IAa6b_saved_parameters_weighed.csv")




#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 10. CREATE MAP WITH UNCERTAIN ECOSYSTEMS ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# – Load Data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
r <- getData("worldclim",var="bio",res=10)
#Bio 1 and Bio12 are mean anual temperature and anual precipitation:
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")


## – Set parameters ####
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Visual Parameters
theme_set(theme_bw()) # set ggplot layout

## – Correct data ####
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Correct for the scaling factor that WorldClim uses
r$Temp = r$Temp/10

# Ecosystems uncertain were mapped for all pixels with MAP > 7.143 MAT + 286 and MAP < –1.469 MAT2 + 81.665 MAT + 475
Con1 = r$Prec > 7.143 * r$Temp + 286
Con2 = r$Prec < -1.469 * r$Temp^2 + 81.665 * r$Temp + 475

plot(Con1)
plot(Con2)
r.combined = Con1 + Con2
r.combined = r.combined >= 2 # select areas that meet both conditions
plot(r.combined) # Uncertain Areas

## Reproject Uncertain Areas 
r.uncertain.repr = projectRaster(from = r.combined, to = raster.example, crs = crs(raster.example), res = res(raster.example)) 
r.uncertain.repr = r.uncertain.repr > 0 # include uncertain areas
plot(r.uncertain.repr)

# Add variable to indicate which cells are uncertain
v.uncertain.repr = as.numeric(getValues(r.uncertain.repr))
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
                #    ggtitle("Current VS present-natural") + # add title
                labs(x="", y="", fill = "Dissimilarity") +
                theme(plot.title = element_text("sans", "plain", "Black", 18, hjust = 0.5),
                      legend.position = "bottom",
                      axis.ticks = element_blank(),
                      axis.text = element_blank(),
                      legend.background = element_blank(),
                      legend.title=element_blank(),
                      legend.text=element_text(size=16)) +
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
ggarrange(p.functional.change, p.veg.change, nrow = 2, labels = "AUTO")
dev.off()

## – End Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
end_time <- Sys.time()
end_time - start_time
