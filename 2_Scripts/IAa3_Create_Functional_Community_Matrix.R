######################################################################!
#                                                                    #
#      IAa5: MAKE FUNCTIONAL GROUP COMMUNITY MATRIX                  #  
#                                                                    #
######################################################################!
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # METADATA
    # This takes functional groups and community matrix and makes a community matrix of functional groups
    # Output: a herbivores traits file that also includes which cluster they belong to
    # time: 18s to run.
    # Status: Ready (2025/02/03)
    # Note: 
    #       We also update the species community matrix so it includes spatial metadata.
    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 1. PREPARING THE DATA ##### 
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####


# – Clear Memory ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
rm(list = ls())  # clean memory (= remove)
graphics.off()  # close graphic windows


# – Setting working directory ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
getwd()
#setwd("/Users/au572919/Research/I_MegaPast2Future/IA_Herbivore_Functional_Groups/1_files/2_component")

# – Setting up libraries ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
library("rworldmap") # visualisation data
library(raster) # GIS
library(sp) # GIS
library(rgeos) # GIS
library(rgdal) # GIS
library(ggplot2)

## – Start Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
start_time <- Sys.time() # Save Time

# – Loading data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
herbivores.4 = read.csv("./3_Output/IAa2_Manual_Functional_Groups/IAa2_herbivores.csv")

community.matrix.current = read.csv("3_Output/IAa1_Create_Community_Matrices/IAa1_current_range_community_matrix.csv")
community.matrix.current = community.matrix.current[,!(colnames(community.matrix.current) %in% c("X","CellNr"))] # remove two redundant variables
community.matrix.current.species = community.matrix.current

community.matrix.present.natural = read.csv("3_Output/IAa1_Create_Community_Matrices/IAa1_present_natural_community_matrix.csv")
community.matrix.present.natural = community.matrix.present.natural[,!(colnames(community.matrix.present.natural) %in% c("X","CellNr"))] # remove two redundant variables
community.matrix.present.natural.species = community.matrix.present.natural

community.matrix.Hempson = read.csv("3_Output/IAa1_Create_Community_Matrices/IAa1_Hempson_community_matrix.csv")
community.matrix.Hempson = community.matrix.Hempson[,!(colnames(community.matrix.Hempson) %in% c("X","CellNr"))] # remove two redundant variables

# Spatial
raster.example = raster("./1_Input/PHYLACINE_V1_2_1/Ranges/Present_natural/Stegodon_orientalis.tif") # example raster from PHYLACINE
biogeo.realms = readOGR("./1_Input/Biogeographical_Realms/udvardy/udvardy_py.shp") # Shapefile of the biogeographical relams https://data-gis.unep-wcmc.org/portal/home/item.html?id=7f3a055f36104f36b6dd7834ebe2cf45

# – Correct data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
herbivores.4$Functional.Group = herbivores.4$Functional.Group.Manual # this line changes the code so we use the manual classifications
functional.group.nr = length(unique(herbivores.4$Functional.Group))
unique(herbivores.4$Functional.Group)

# – Parameters ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
theme_set(theme_bw()) # set ggplot layout


## – Directories  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
if(!dir.exists(paste0("./3_Output/IAa3_Create_Functional_Community_Matrix/"))){
  dir.create(paste0("./3_Output/IAa3_Create_Functional_Community_Matrix/"))}



#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 2. MAKE THE CURRENT FUNCTIONAL GROUP COMMUNITY MATRIX ##### 
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# – Current ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

functional.community.matrix.current = matrix(NA, nrow = nrow(community.matrix.current), ncol = length(unique(herbivores.4$Functional.Group)) ) # create matrix

for (i.group in 1:length(unique(herbivores.4$Functional.Group))) {
  temp = which(colnames(community.matrix.current) %in% herbivores.4[which(herbivores.4$Functional.Group == unique(herbivores.4$Functional.Group)[i.group]),"Binomial"])
  functional.community.matrix.current[,i.group] = apply(data.frame(community.matrix.current[,temp]), 1, sum)
}

colnames(functional.community.matrix.current) = paste("Functional.Group", unique(herbivores.4$Functional.Group), sep =".")
str(functional.community.matrix.current)
range(functional.community.matrix.current)

# – Present-natural ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

functional.community.matrix.present.natural = matrix(NA, nrow = nrow(community.matrix.present.natural), ncol = length(unique(herbivores.4$Functional.Group)))  # create matrix

for (i.group in 1:length(unique(herbivores.4$Functional.Group))) {
  temp = which(colnames(community.matrix.present.natural) %in% herbivores.4[which(herbivores.4$Functional.Group == unique(herbivores.4$Functional.Group)[i.group]),"Binomial"])
  functional.community.matrix.present.natural[,i.group] = apply(data.frame(community.matrix.present.natural[,temp]), 1, sum)
}

colnames(functional.community.matrix.present.natural) = paste("Functional.Group", unique(herbivores.4$Functional.Group), sep =".")
str(functional.community.matrix.present.natural)
range(functional.community.matrix.present.natural)

# – Hempson ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

functional.community.matrix.Hempson = matrix(NA, nrow = nrow(community.matrix.Hempson), ncol = length(unique(herbivores.4$Functional.Group)))  # create matrix

for (i.group in 1:length(unique(herbivores.4$Functional.Group))) {
  temp = which(colnames(community.matrix.Hempson) %in% herbivores.4[which(herbivores.4$Functional.Group == unique(herbivores.4$Functional.Group)[i.group]),"Binomial"])
  functional.community.matrix.Hempson[,i.group] = apply(data.frame(community.matrix.Hempson[,temp]), 1, sum)
}

colnames(functional.community.matrix.Hempson) = paste("Functional.Group", unique(herbivores.4$Functional.Group), sep =".")
str(functional.community.matrix.Hempson)
range(functional.community.matrix.Hempson)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 3. MAKE THE COMBINED GROUP COMMUNITY MATRIX ##### 
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

      # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #  
      #   Here, we combine the current and the present-natural community matrices into one big matrix. 
      #   We do so, to enable later analyses where we will calculate the distance between current and present-natural cells,
      #   As well as the distances within these scenarios.
      #
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# – FUNCTIONAL GROUPS ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

functional.community.matrix.combined = data.frame(rbind(functional.community.matrix.current, functional.community.matrix.present.natural))
str(functional.community.matrix.combined)

# turn community matrix into df to add new columwn
functional.community.matrix.current = data.frame(functional.community.matrix.current)
functional.community.matrix.present.natural = data.frame(functional.community.matrix.present.natural)
functional.community.matrix.Hempson = data.frame(functional.community.matrix.Hempson)

functional.community.matrix.current$Scenario = "Current"
functional.community.matrix.present.natural$Scenario = "Present-Natural"
functional.community.matrix.Hempson$Scenario = "Hempson"

functional.community.matrix.current$Original.Cell = rownames(functional.community.matrix.current)
functional.community.matrix.present.natural$Original.Cell = rownames(functional.community.matrix.present.natural)
functional.community.matrix.Hempson$Original.Cell = rownames(functional.community.matrix.Hempson)

# Add columns to combined data frame.
meta.variables = data.frame(rbind(functional.community.matrix.current[,c("Scenario", "Original.Cell")],functional.community.matrix.present.natural[,c("Scenario", "Original.Cell")]))
functional.community.matrix.combined = cbind(functional.community.matrix.combined,meta.variables)
str(functional.community.matrix.combined)

# change data type
functional.community.matrix.combined$Scenario = as.factor(functional.community.matrix.combined$Scenario)
functional.community.matrix.combined$Original.Cell = as.factor(functional.community.matrix.combined$Original.Cell)
functional.community.matrix.Hempson$Original.Cell = as.factor(functional.community.matrix.Hempson$Original.Cell)

#check data
colnames(functional.community.matrix.combined)
dim(functional.community.matrix.combined)
str(functional.community.matrix.combined)

colnames(functional.community.matrix.Hempson)
dim(functional.community.matrix.Hempson)
str(functional.community.matrix.Hempson)

# – SPECIES ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
community.matrix.combined.species = data.frame(rbind(community.matrix.current.species, community.matrix.present.natural.species))
#str(community.matrix.combined.species)

# turn community matrix into df to add new columwn
community.matrix.current.species = data.frame(community.matrix.current.species)
community.matrix.present.natural.species = data.frame(community.matrix.present.natural.species)

community.matrix.current.species$Scenario = "Current"
community.matrix.present.natural.species$Scenario = "Present-Natural"

community.matrix.current.species$Original.Cell = rownames(community.matrix.current.species)
community.matrix.present.natural.species$Original.Cell = rownames(community.matrix.present.natural.species)

# Add columns to combined data frame.
meta.variables = data.frame(rbind(community.matrix.current.species[,c("Scenario", "Original.Cell")],community.matrix.present.natural.species[,c("Scenario", "Original.Cell")]))
community.matrix.combined.species = cbind(community.matrix.combined.species,meta.variables)
str(community.matrix.combined.species)

# change data type
community.matrix.combined.species$Scenario = as.factor(community.matrix.combined.species$Scenario)
community.matrix.combined.species$Original.Cell = as.factor(community.matrix.combined.species$Original.Cell)

#check data
colnames(community.matrix.combined.species)
dim(community.matrix.combined.species)
str(community.matrix.combined.species)

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 4. ADD SPATIAL METADATA ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# – Make Countries Polygon ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
data("countriesCoarse")
countries = (countriesCoarse)
countries@data$continent[which(countries@data$NAME == "French Guiana")] = "South America"
countries = spTransform(countries, CRSobj = crs(raster.example))
countries.poly = (countries[(countries@data$continent %in% c("Eurasia", "Africa", "Australia", "North America", "South America")),])

countries.raster = rasterize(countries.poly,raster.example, field = 1)
#plot(countries.raster)
#plot(countries.poly, add = T)

# – Dissolve Biogeographical realm ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
biogeo.realms = spTransform(biogeo.realms, proj4string(countries.raster))
biogeo.realms.dissolved <- gUnaryUnion(biogeo.realms, id = biogeo.realms@data$realmnum)
plot(biogeo.realms.dissolved)
biogeo.realms.raster = rasterize(biogeo.realms.dissolved,raster.example, background = NA)
biogeo.realms.raster = mask(biogeo.realms.raster,countries.raster)
#plot(biogeo.realms.raster)

# select all cell (at continent borders) that contain herbivores but that do not belong to a biogeographical realm and set them to NA (they are noise)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # NOTE: 
      #   Because of a mismatch between layers we had various cells at the edge of continents that did contain herbivores but that did not belong to a biogeographical realm.
      #   As a solution we clipped everything to the countries.raster file, so the cells at the edge would be removed.
      #   This removed most of the cells but not all.Instead, we just manually set all these borders cells to NA (after visually checking on map)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ––– Script to check for problematic cells ––––
#
# temp = df.master[which(is.na(df.master$Biogeo.Realm) & df.master$Presence == 1),]
# ggplot() +
#   geom_raster(data = temp, aes(x=Longitude, y=Latitude, fill="red")) +
#   geom_polygon(data = countries.poly, 
#                aes(x = long, 
#                    y = lat, 
#                    group = group),
#                fill = NA, 
#                colour = 'gray10', 
#                size = 0.1) +
#   theme_bw() +
#   coord_fixed() + 
#   ggtitle("Present-Natural") + # add title
#   theme(plot.title = element_text(hjust = 0.5)) # centre title 

# Biogeographical Realms 
functional.community.matrix.combined$Biogeo.realm = c(getValues(biogeo.realms.raster), getValues(biogeo.realms.raster)) # add Biogeo.realm variable
functional.community.matrix.combined$Countries.raster = c(getValues(countries.raster), getValues(countries.raster)) # add countries.raster variable
functional.community.matrix.combined[which(is.na(functional.community.matrix.combined$Biogeo.realm) | (!(functional.community.matrix.combined$Countries.raster == 1))),] = NA # Turn cells that do not occur in a biogeographical realm or in a the countries.raster mask to NA.

functional.community.matrix.Hempson$Biogeo.realm = c(getValues(biogeo.realms.raster)) # add Biogeo.realm variable
functional.community.matrix.Hempson$Countries.raster = c(getValues(countries.raster)) # add countries.raster variable
functional.community.matrix.Hempson[which(is.na(functional.community.matrix.Hempson$Biogeo.realm) | (!(functional.community.matrix.Hempson$Countries.raster == 1))),] = NA # Turn cells that do not occur in a biogeographical realm or in a the countries.raster mask to NA.

community.matrix.combined.species$Biogeo.realm = c(getValues(biogeo.realms.raster), getValues(biogeo.realms.raster)) # add Biogeo.realm variable
community.matrix.combined.species$Countries.raster = c(getValues(countries.raster), getValues(countries.raster)) # add countries.raster variable
community.matrix.combined.species[which(is.na(community.matrix.combined.species$Biogeo.realm) | (!(community.matrix.combined.species$Countries.raster == 1))),] = NA # Turn cells that do not occur in a biogeographical realm or in a the countries.raster mask to NA.

# – Name Biogeographical realms ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

functional.community.matrix.combined$Biogeo.realm = as.character(functional.community.matrix.combined$Biogeo.realm) # turn biogeo realms in character string to add names
functional.community.matrix.combined[which(functional.community.matrix.combined$Biogeo.realm == 2),"Biogeo.realm"] = "Nearctic"
functional.community.matrix.combined[which(functional.community.matrix.combined$Biogeo.realm == 3),"Biogeo.realm"] = "Palaearctic"
functional.community.matrix.combined[which(functional.community.matrix.combined$Biogeo.realm == 4),"Biogeo.realm"] = "Africotropical"
functional.community.matrix.combined[which(functional.community.matrix.combined$Biogeo.realm == 5),"Biogeo.realm"] = "Indomalayan"
functional.community.matrix.combined[which(functional.community.matrix.combined$Biogeo.realm == 6),"Biogeo.realm"] = "Australasian"
functional.community.matrix.combined[which(functional.community.matrix.combined$Biogeo.realm == 7),"Biogeo.realm"] = "Australasian"
functional.community.matrix.combined[which(functional.community.matrix.combined$Biogeo.realm == 8),"Biogeo.realm"] = "Oceanian"
functional.community.matrix.combined[which(functional.community.matrix.combined$Biogeo.realm == 9),"Biogeo.realm"] = "Neotropical"
functional.community.matrix.combined$Biogeo.realm = as.factor(functional.community.matrix.combined$Biogeo.realm) # turn names into factors

functional.community.matrix.Hempson$Biogeo.realm = as.character(functional.community.matrix.Hempson$Biogeo.realm) # turn biogeo realms in character string to add names
functional.community.matrix.Hempson[which(functional.community.matrix.Hempson$Biogeo.realm == 2),"Biogeo.realm"] = "Nearctic"
functional.community.matrix.Hempson[which(functional.community.matrix.Hempson$Biogeo.realm == 3),"Biogeo.realm"] = "Palaearctic"
functional.community.matrix.Hempson[which(functional.community.matrix.Hempson$Biogeo.realm == 4),"Biogeo.realm"] = "Africotropical"
functional.community.matrix.Hempson[which(functional.community.matrix.Hempson$Biogeo.realm == 5),"Biogeo.realm"] = "Indomalayan"
functional.community.matrix.Hempson[which(functional.community.matrix.Hempson$Biogeo.realm == 6),"Biogeo.realm"] = "Australasian"
functional.community.matrix.Hempson[which(functional.community.matrix.Hempson$Biogeo.realm == 7),"Biogeo.realm"] = "Australasian"
functional.community.matrix.Hempson[which(functional.community.matrix.Hempson$Biogeo.realm == 8),"Biogeo.realm"] = "Oceanian"
functional.community.matrix.Hempson[which(functional.community.matrix.Hempson$Biogeo.realm == 9),"Biogeo.realm"] = "Neotropical"
functional.community.matrix.Hempson$Biogeo.realm = as.factor(functional.community.matrix.Hempson$Biogeo.realm) # turn names into factors

community.matrix.combined.species$Biogeo.realm = as.character(community.matrix.combined.species$Biogeo.realm) # turn biogeo realms in character string to add names
community.matrix.combined.species[which(community.matrix.combined.species$Biogeo.realm == 2),"Biogeo.realm"] = "Nearctic"
community.matrix.combined.species[which(community.matrix.combined.species$Biogeo.realm == 3),"Biogeo.realm"] = "Palaearctic"
community.matrix.combined.species[which(community.matrix.combined.species$Biogeo.realm == 4),"Biogeo.realm"] = "Africotropical"
community.matrix.combined.species[which(community.matrix.combined.species$Biogeo.realm == 5),"Biogeo.realm"] = "Indomalayan"
community.matrix.combined.species[which(community.matrix.combined.species$Biogeo.realm == 6),"Biogeo.realm"] = "Australasian"
community.matrix.combined.species[which(community.matrix.combined.species$Biogeo.realm == 7),"Biogeo.realm"] = "Australasian"
community.matrix.combined.species[which(community.matrix.combined.species$Biogeo.realm == 8),"Biogeo.realm"] = "Oceanian"
community.matrix.combined.species[which(community.matrix.combined.species$Biogeo.realm == 9),"Biogeo.realm"] = "Neotropical"
community.matrix.combined.species$Biogeo.realm = as.factor(community.matrix.combined.species$Biogeo.realm) # turn names into factors

#check data
colnames(functional.community.matrix.combined)
dim(functional.community.matrix.combined)
str(functional.community.matrix.combined)

colnames(functional.community.matrix.Hempson)
dim(functional.community.matrix.Hempson)
str(functional.community.matrix.Hempson)

colnames(community.matrix.combined.species)
dim(community.matrix.combined.species)
str(community.matrix.combined.species)

# – Add Coordinates ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

## Current and PN
df.master <- rasterToPoints(raster.example) # create points with spatial information
df.master[,3] = getValues(countries.raster) # add variable showing which cells fall within the countriesraster folder.
df.master = rbind(df.master,df.master) # duplicate the file for current and present-natural points
colnames(df.master) <- c("Longitude", "Latitude", "countries.raster") # Make appropriate column headings
df.master = data.frame(df.master)

functional.community.matrix.combined = cbind(functional.community.matrix.combined, df.master[,c("Longitude", "Latitude")]) # bind df.master to the combined community matrix
community.matrix.combined.species = cbind(community.matrix.combined.species, df.master[,c("Longitude", "Latitude")]) # bind df.master to the combined community matrix

## Hempson
functional.community.matrix.Hempson = data.frame(functional.community.matrix.Hempson)
functional.community.matrix.Hempson$Original.Cell = rownames(functional.community.matrix.Hempson) # Add Cell NR

df.master.2 <- rasterToPoints(raster.example) # create points with spatial information
df.master.2[,3] = getValues(countries.raster) # add variable showing which cells fall within the countriesraster folder.
colnames(df.master.2) <- c("Longitude", "Latitude", "countries.raster") #Make appropriate column headings
df.master.2 = data.frame(df.master.2)

functional.community.matrix.Hempson = cbind(functional.community.matrix.Hempson, df.master.2[,c("Longitude", "Latitude")])

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 5. CREATE DIFFERENT COMMUNITY MATRICES ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# Select the functional group columns
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
columns.with.functional.groups = which(grepl("Functional.Group.", colnames(functional.community.matrix.combined)))

# – Select cells with herbivores ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
cells.with.herbivores = which(apply(functional.community.matrix.combined[,columns.with.functional.groups],1,sum) > 0) # select cells where herbivores occur
functional.community.matrix.combined$Herbivores.Present = 0 # Add this as a variable to the df
functional.community.matrix.combined[cells.with.herbivores,"Herbivores.Present"] = 1

cells.with.herbivores.Hempson = which(apply(functional.community.matrix.Hempson[,columns.with.functional.groups],1,sum) > 0) # select cells where herbivores occur
functional.community.matrix.Hempson$Herbivores.Present = 0 # Add this as a variable to the df
functional.community.matrix.Hempson[cells.with.herbivores.Hempson,"Herbivores.Present"] = 1

columns.with.species = which(grepl("_", colnames(community.matrix.combined.species))) # select all colnames with an underscore (only species names have these)
cells.with.herbivores = which(apply(community.matrix.combined.species[,columns.with.species],1,sum) > 0) # select cells where herbivores occur
community.matrix.combined.species$Herbivores.Present = 0 # Add this as a variable to the df
community.matrix.combined.species[cells.with.herbivores,"Herbivores.Present"] = 1

# – Relative abundance community matrix ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

      # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #  
      #   This community matrix shows the proportion of the assemblages that is occupied by the functional type.
      #   The proportions can range from 0 (not present) to 1 (only functional type present).
      #
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

functional.community.matrix.relative = functional.community.matrix.combined # make a copy of functional.community.matrix
functional.community.matrix.relative[cells.with.herbivores,columns.with.functional.groups] = functional.community.matrix.relative[cells.with.herbivores,columns.with.functional.groups] / apply(functional.community.matrix.relative[cells.with.herbivores,columns.with.functional.groups], 1, sum) # calculate relative abundance values expressed as a percentage  

functional.community.matrix.relative.Hempson = functional.community.matrix.Hempson # make a copy of functional.community.matrix
functional.community.matrix.relative.Hempson[cells.with.herbivores.Hempson,columns.with.functional.groups] = functional.community.matrix.relative.Hempson[cells.with.herbivores.Hempson,columns.with.functional.groups] / apply(functional.community.matrix.relative.Hempson[cells.with.herbivores.Hempson,columns.with.functional.groups], 1, sum) # calculate relative abundance values expressed as a percentage  

# check data
colnames(functional.community.matrix.relative)
range(functional.community.matrix.relative[,columns.with.functional.groups], na.rm = T)
dim(functional.community.matrix.relative)
str(functional.community.matrix.relative)

colnames(functional.community.matrix.relative.Hempson)
range(functional.community.matrix.relative.Hempson[,columns.with.functional.groups], na.rm = T)
dim(functional.community.matrix.relative.Hempson)
str(functional.community.matrix.relative.Hempson)

# – Presence/Absence Community Matrix ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

      # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #  
      #   This community matrix shows whether a functional type is present (1) or absent (0) in an assemblage.
      #
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

functional.community.matrix.presence.absence = functional.community.matrix.combined
functional.community.matrix.presence.absence = functional.community.matrix.presence.absence[which(functional.community.matrix.presence.absence$Herbivores.Present == 1),]

temp = functional.community.matrix.presence.absence[,columns.with.functional.groups]
temp[temp > 0] = 1 

functional.community.matrix.presence.absence[,columns.with.functional.groups] = temp

# check data
colnames(functional.community.matrix.presence.absence)
range(functional.community.matrix.presence.absence[,1:10])
dim(functional.community.matrix.presence.absence)

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 5. CREATE WEIGHED COMMUNITY MATRIX ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####


    # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  
    #   Adding weights in the Analog Distance function does not work when using a Squared Chord Distance. 
    #   Therefore we are adding weights manually. 
    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create Weights
columns.with.functional.groups = which(grepl("Functional.Group.", colnames(functional.community.matrix.relative)))
v.weights = rep(NA,length(columns.with.functional.groups)) ## create matrix

# v.weights[which(grepl("meso", colnames(functional.community.matrix.relative)))] = 1 # mesoherbivore weight
# v.weights[which(grepl("macro", colnames(functional.community.matrix.relative)))] = 1.5 # macroherbivore weight
# v.weights[which(grepl("mega", colnames(functional.community.matrix.relative)))] = 2 # megaherbivore weight


matches.meso = grep(x = herbivores.4$Functional.Group , pattern = "meso")
weight.meso = (mean(herbivores.4[matches.meso,"Mass.g"])) 

matches.macro = grep(x = herbivores.4$Functional.Group , pattern = "macro")
weight.macro = (mean(herbivores.4[matches.macro,"Mass.g"])) 

matches.mega = grep(x = herbivores.4$Functional.Group , pattern = "mega")
weight.mega = (mean(herbivores.4[matches.mega,"Mass.g"])) 


v.weights[which(grepl("meso", colnames(functional.community.matrix.relative)))] = log10(weight.meso)/log10(weight.meso) # mesoherbivore weight
v.weights[which(grepl("macro", colnames(functional.community.matrix.relative)))] = log10(weight.macro)/log10(weight.meso) # macroherbivore weight
v.weights[which(grepl("mega", colnames(functional.community.matrix.relative)))] = log10(weight.mega)/log10(weight.meso) # megaherbivore weight

## weight the columns according to weights
df.temp =  data.frame(as.matrix(functional.community.matrix.relative[,columns.with.functional.groups]) %*% diag(v.weights)) 
df.temp.Hempson =  data.frame(as.matrix(functional.community.matrix.relative.Hempson[,columns.with.functional.groups]) %*% diag(v.weights)) 


## Restandardise the values so they become percentages again.
df.temp = t(apply(df.temp, 1, function(x){ x/sum(x, na.rm = T)}))
dim(df.temp)
range(df.temp, na.rm = T)

df.temp.Hempson = t(apply(df.temp.Hempson, 1, function(x){ x/sum(x, na.rm = T)}))
dim(df.temp.Hempson)
range(df.temp.Hempson, na.rm = T)

## Overwrite values
functional.community.matrix.relative.weighed = functional.community.matrix.relative
functional.community.matrix.relative.weighed[,columns.with.functional.groups] = df.temp
str(functional.community.matrix.relative.weighed)
rm(df.temp) # remove clutter

functional.community.matrix.relative.Hempson.weighed = functional.community.matrix.relative.Hempson
functional.community.matrix.relative.Hempson.weighed[,columns.with.functional.groups] = df.temp.Hempson
str(functional.community.matrix.relative.Hempson.weighed)
rm(df.temp.Hempson) # remove clutter

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 6. SAVE DATA ##### 
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

write.csv(functional.community.matrix.combined, "./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_abundances.csv", row.names = F)
write.csv(functional.community.matrix.presence.absence, "./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_presence_absence.csv", row.names = F)
write.csv(functional.community.matrix.relative, "./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_relative_abundances.csv", row.names = F)
write.csv(functional.community.matrix.relative.weighed, "./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_relative_abundances_weighed.csv", row.names = F)
write.csv(functional.community.matrix.Hempson, "./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_Hempson_abundances.csv", row.names = F)
write.csv(functional.community.matrix.relative.Hempson, "./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_Hempson_relative_abundances.csv", row.names = F)
write.csv(functional.community.matrix.relative.Hempson.weighed, "./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_Hempson_relative_abundances_weighed.csv", row.names = F)
write.csv(community.matrix.combined.species, "./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_community_matrix_combined_species.csv", row.names = F)


## – End Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
end_time <- Sys.time()
end_time - start_time

