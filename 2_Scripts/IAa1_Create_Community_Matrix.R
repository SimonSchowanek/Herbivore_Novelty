#######################################################################!
#                                                                    #
#      1Aa1 CREATING THE COMMUNITY MATRIX                            #                                            
#                                                                    #
#####################################################################!

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # METADATA
    # This file creates the community matrices required in the later analysis
    #
    # Output: one current and one present natural community matrix
    # Time: needs 2 min to run.
    # Status: Ready (2025/02/03)
    # NOTES: Dasypus bellis has been removed based on diet.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 1. PREPARING THE DATA ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# – Clear Memory ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
rm(list = ls())  # clean memory (= remove)
graphics.off()  # close graphic windows


# – Setting working directory ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
getwd()
#setwd("/Users/au572919/Research/I_MegaPast2Future/IA_Herbivore_Functional_Groups/1_files/2_Component")

# – Set up libraries ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
library(sp)
library(rgdal)
library("raster")

## – Start Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
start_time <- Sys.time() # Save Time


# – Load data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
herbivores = read.csv("/Users/au572919/Research/I_MegaPast2Future/HerbiTraits/Data/Version_1.1/HerbiTraits_1.1/HerbiTraits_1.1.csv") # herbitraits (https://github.com/MegaPast2Future/HerbiTraits)
raster_base = raster("/Users/au572919/Datasets/Phylacine.1.2/Data/Ranges/Current/Abditomys_latidens.tif") # example raster from PHYLACINE (https://megapast2future.github.io/PHYLACINE_1.2/)
spatial.metadata = read.csv("/Users/au572919/Datasets/Phylacine.1.2/Data/Ranges/Spatial_metadata.csv") # spatial metadata from PHYLACINE

# – Correct data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
herbivores = herbivores[which(herbivores$Class == "Mammalia"),] # remove birds
herbivores = herbivores[which(herbivores$Domesticated == 0),] # remove domesticated herbivores
herbivores = herbivores[order(herbivores$Binomial),]

herbivores$Binomial = gsub(" ", "_", herbivores$Binomial)



## – Directories  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
if(!dir.exists(paste0("./3_Output"))){
  dir.create(paste0("./3_Output"))}


if(!dir.exists(paste0("./3_Output/IAa1_Create_Community_Matrices/"))){
  dir.create(paste0("./3_Output/IAa1_Create_Community_Matrices/"))}


#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 2. CREATE LIST OF UNIQUE SPECIES ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# remove Dassipus bellus ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# herbivores = herbivores[which(!(herbivores$Binomial == "Dasypus_bellus")),]

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # NOTE: 
      #   We remove Dasypus_bellus because it should not be in the list (mostly omnivorous)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      

# create extant and extinct list ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
herbivores.extinct = herbivores[which(herbivores$Species.Status %in% c("Extinct before 1500 CE","Extinct after 1500 CE")),]
herbivores.extant = herbivores[which(!(herbivores$Species.Status %in% c("Extinct before 1500 CE","Extinct after 1500 CE"))),]

# remove species that are not mapped by IUCN ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
species.with.ranges = spatial.metadata[which(!((spatial.metadata$Number.Cells.Current.Range == 0) & (spatial.metadata$Number.Cells.Present.Natural.Range == 0))),"Binomial.1.2"]
length(species.with.ranges)

herbivores.2 = herbivores[which(herbivores$Binomial %in% species.with.ranges),] # check which species occur in species.with.ranges

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # NOTE: 
        #   We remove Piliocolobus_pennantii & Gazella Marica as their ranges have not been mapped by IUCN
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a list with unique species
List_Of_Unique_Species = unique(paste(herbivores.2$Binomial, ".tif", sep = "")) # species list
List_Of_Unique_Species = as.vector(List_Of_Unique_Species)
List_Of_Unique_Species = List_Of_Unique_Species[which(!(List_Of_Unique_Species == ""))]

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 3. CREATE PRESENT NATURAL COMMUNITY MATRIX ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

herbivores.2[which(!(herbivores.2$Present_Natural_Map_Name %in% List_Of_Unique_Species)),] # identify how many species lack a present natural range


# preparing for the loop
directory = "/Users/au572919/Datasets/Phylacine.1.2/Data/Ranges/Present_natural/" # creating directory for intput 

Nr_Of_Runs = length(List_Of_Unique_Species)
Community_Matrix_PN = 1:ncell(raster_base)

Map_Name = paste(herbivores.2$Binomial, ".tif", sep = "")

# Making the community matrix
for(i in 1:Nr_Of_Runs)
{
  #defining the rastername
  binomial = List_Of_Unique_Species[i]
  
  if (binomial %in% Map_Name) {
    
    raster_name = paste(directory,binomial,sep ="")    
    r = raster(raster_name) #reading the data
    vector = getValues(r) #making a vector
    Community_Matrix_PN = cbind(Community_Matrix_PN, vector) #pasting it to the existing vector
    
  } else {
    
    absent_species = rep(x = 0, ncell(raster_base))
    Community_Matrix_PN = cbind(Community_Matrix_PN, absent_species) #pasting it to the existing vector
    
  }
  
  #update on progress  
  progress = (i / Nr_Of_Runs) * 100
  print(progress)
  
}

colnames(Community_Matrix_PN)[1] = "CellNr" # define column names
colnames(Community_Matrix_PN)[2: (Nr_Of_Runs + 1)] = as.character(herbivores.2$Binomial) # define column names

write.csv( x= Community_Matrix_PN,
           file = "3_Output/IAa1_Create_Community_Matrices/IAa1_present_natural_community_matrix.csv",
           row.names = FALSE)

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 4. CREATE CURRENT RANGES COMMUNITY MATRIX ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# identify how many species lack a present natural range
herbivores.2[which(!(herbivores.2$Present_Natural_Map_Name %in% List_Of_Unique_Species)),]

# preparing for the loop
directory = "/Users/au572919/Datasets/Phylacine.1.2/Data/Ranges/Current/" # creating directory for output 
Nr_Of_Runs = length(List_Of_Unique_Species)
Community_Matrix_Current = 1:ncell(raster_base)

Map_Name = paste(herbivores.2$Binomial, ".tif", sep = "")

# Making the community matrix
for(i in 1:Nr_Of_Runs)
{
  #defining the rastername
  binomial = List_Of_Unique_Species[i]
  
  if (binomial %in% Map_Name) {
    
    raster_name = paste(directory,binomial,sep ="")    
    r = raster(raster_name) #reading the data
    vector = getValues(r) #making a vector
    Community_Matrix_Current = cbind(Community_Matrix_Current, vector) #pasting it to the existing vector
    
  } else {
    
    absent_species = rep(x = 0, ncell(raster_base))
    Community_Matrix_Current = cbind(Community_Matrix_Current, absent_species) #pasting it to the existing vector
    
  }
  
  #update on progress  
  progress = (i / Nr_Of_Runs) * 100
  print(progress)
  
}

colnames(Community_Matrix_Current)[1] = "CellNr" # define column names
colnames(Community_Matrix_Current)[2: (Nr_Of_Runs + 1)] = as.character(herbivores.2$Binomial) # define column names

write.csv( x= Community_Matrix_Current,
           file = "3_Output/IAa1_Create_Community_Matrices/IAa1_current_range_community_matrix.csv",
           row.names = F)

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 5. HEMPSON COMMUNITY MATRIX ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

      # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #  
      #   We caclulate how different the Hempson et al. (2015) herbivomes are from each other. Therefore we need a community matrix that 
      #   is comparable to the community matrices that Hempson made. Hempson et al. (2015) use the current ranges for African herbivores. However,
      #   they modify the range of Giraffa_camelopardalis, Diceros_bicornis, Ceratotherium_simum and Equus_zebra so they better reflect the historical distributions.
      #   We (visually) compared the current and present-natural ranges of these species, and the Hempson modified range best match the present-natural ranges.
      #   Our solution is to take the current ranges of all species, but to replace the ranges of the four aforementioned species with their present-natural ranges.
      #
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Community_Matrix_Hempson = Community_Matrix_Current

## Overwrite those four species
Community_Matrix_Hempson[,which(colnames(Community_Matrix_Hempson) == "Giraffa_camelopardalis")] = Community_Matrix_PN[,which(colnames(Community_Matrix_PN) == "Giraffa_camelopardalis")]
Community_Matrix_Hempson[,which(colnames(Community_Matrix_Hempson) == "Diceros_bicornis")] = Community_Matrix_PN[,which(colnames(Community_Matrix_PN) == "Diceros_bicornis")]
Community_Matrix_Hempson[,which(colnames(Community_Matrix_Hempson) == "Ceratotherium_simum")] = Community_Matrix_PN[,which(colnames(Community_Matrix_PN) == "Ceratotherium_simum")]
Community_Matrix_Hempson[,which(colnames(Community_Matrix_Hempson) == "Equus_zebra")] = Community_Matrix_PN[,which(colnames(Community_Matrix_PN) == "Equus_zebra")]


max(abs((Community_Matrix_Hempson) - Community_Matrix_Current)) # check if there are differences between the matrices

write.csv( x= Community_Matrix_Hempson,
           file = "3_Output/IAa1_Create_Community_Matrices/IAa1_Hempson_community_matrix.csv",
           row.names = F)

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 6. SAVE ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

write.csv( x= herbivores.2,
           file = "3_Output/IAa1_Create_Community_Matrices/IAa1_herbivores_2.csv", row.names = F)

write.csv( x= herbivores.extinct,
           file = "3_Output/IAa1_Create_Community_Matrices/IAa1_herbivores_extinct.csv", row.names = F)

write.csv( x= herbivores.extant,
           file = "3_Output/IAa1_Create_Community_Matrices/IAa1_herbivores_extant.csv", row.names = F)

## – End Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
end_time <- Sys.time()
end_time - start_time
