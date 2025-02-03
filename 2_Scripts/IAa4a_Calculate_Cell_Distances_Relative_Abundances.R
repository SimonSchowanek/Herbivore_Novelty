######################################################################!
#                                                                    #
#      IAa4: MAKE CELL DISTANCE MATRIX BASED ON RELATIVE ABUNDANCES  #  
#                                                                    #
######################################################################!

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # METADATA
    # This script calculates the community composition distance between raster cells.
    #
    # Output: a distance object and a distance matrix
    # time: 10 min
    # Status: Ready (2025/02/03)
    # Note: 
    #      
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

# – Setting up libraries ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
library(dplyr) # for data cleaningx
library(analogue) # to perform the distance calculation
library(bigmemory) # To work with large datasets within R


## – Start Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
start_time <- Sys.time() # Save Time


# – Load Data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
functional.community.matrix = read.csv("3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_relative_abundances.csv")
functional.community.matrix = functional.community.matrix[,!(colnames(functional.community.matrix) %in% c("X","CellNr"))] # remove two redundant variables


# Check the data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

## Combined Data
functional.community.matrix[1:10,1:10]
head(functional.community.matrix)
colnames(functional.community.matrix)
dim(functional.community.matrix)


## – Directories  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
if(!dir.exists(paste0("./3_Output/IAa4_Calculate_Cell_Distances/"))){
  dir.create(paste0("./3_Output/IAa4_Calculate_Cell_Distances/"))}



#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 2. CALCULATE DISTANCES ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# – Select cells with herbivores----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

functional.community.matrix.relative.subsample = functional.community.matrix[which(functional.community.matrix$Herbivores.Present == 1),]
colnames(functional.community.matrix.relative.subsample)
dim(functional.community.matrix.relative.subsample)

# – Select columns with functional groups----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
columns.with.functional.groups = which(grepl("Functional.Group.", colnames(functional.community.matrix.relative.subsample)))


# – Calculate distances between cells ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

      # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #  
      #   Here, we calculate the squared chord distance. The squared chord distance is a distance metric that is commonly
      #   used in palynology to determine the dissimilarity between ecological assemblages. It calculates a dissimilarity
      #   value ranging from zero to two, based on the functional types (or more commonly species) that two assemblages have 
      #   in common. A value of zero means both assemblages share all functional types, in the same proportions (here, relative 
      #   richness). A value of two means that the assemblages have no functional types in common.
      #
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system.time( # 165s
  dist.matrix.abundance.rel.distance.SQchord <- distance(functional.community.matrix.relative.subsample[,columns.with.functional.groups], method = "SQchord", weights = NULL, R = NULL, dist = T, fast = T))

# Check data
str(dist.matrix.abundance.rel.distance.SQchord)

# Save data
system.time( # 125.983s
  save(dist.matrix.abundance.rel.distance.SQchord, file = "./3_Output/IAa4_Calculate_Cell_Distances/IAa4_dist_matrix_abundance_rel_distance_SQchord"))



# – Convert distances to matrix ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

system.time( # 733.537s
  dist.matrix.abundance.rel.distance.SQchord.matrix <- as.matrix(dist.matrix.abundance.rel.distance.SQchord)) # takes less than an hour

# Check data
str(dist.matrix.abundance.rel.distance.SQchord.matrix)

# Remove data
rm(dist.matrix.abundance.rel.distance.SQchord)

# Save data
system.time( # 125.983s
  save(dist.matrix.abundance.rel.distance.SQchord.matrix, file = "./3_Output/IAa4_Calculate_Cell_Distances/IAa4_dist_matrix_abundance_rel_distance_SQchord_matrix"))


system.time( # 3157.486 s
   write.table(dist.matrix.abundance.rel.distance.SQchord.matrix, 
               "./3_Output/IAa4_Calculate_Cell_Distances/IAa4_dist_matrix_abundance_perc_SQchord_matrix.csv",
               row.names = F,
               col.names = F,
               sep = ","))

# – Convert matrix to big matrix ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Load matrix
system.time( #346.015
  dist.matrix.abundance.rel.distance.SQchord.matrix.big <- as.big.matrix(dist.matrix.abundance.rel.distance.SQchord.matrix,
                                                                         type = "double",
                                                                         separated = F,
                                                                         backingfile = "IAa4_big_distance_matrix.bin",
                                                                         descriptorfile = "IAa4_big_distance_matrix.desc",
                                                                         backingpath = "./3_Output/IAa4_Calculate_Cell_Distances/"))

# remove data
rm(dist.matrix.abundance.rel.distance.SQchord.matrix)

# Check data
str(dist.matrix.abundance.rel.distance.SQchord.matrix.big)
dist.matrix.abundance.rel.distance.SQchord.matrix.big[1:10,1:10]
#range(dist.matrix.abundance.rel.distance.SQchord.matrix.big)

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 3. CALCULATE MATRICES FOR CODE TESTING (OPTIONAL) ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

      # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #  
      #   In this sections we save three separate versions of the distance matrix. The large distance matrix is the 
      #   complete version. The medium and the small version as smaller subsets of the large version, based on random samples.
      #   They are included to have an object to test our code with, as calculations using the large distance matrix may take several hours.
      #   The small and medium distances matrices have no importance for further analyses.
      #
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
little.sample = sample(1:nrow(functional.community.matrix.relative.subsample),100, replace = F)
functional.community.matrix.relative.little.subsample = functional.community.matrix.relative.subsample[little.sample,]
# 
# system.time( # 165s
functional.community.matrix.relative.little.subsample.SQchord <- distance(functional.community.matrix.relative.little.subsample[,columns.with.functional.groups], method = "SQchord", weights = NULL, R = NULL, dist = T, fast = T)
save(functional.community.matrix.relative.little.subsample.SQchord,file = "./3_Output/IAa4_Calculate_Cell_Distances/IAa4_dist_matrix_abundance_rel_distance_SQchord_matrix_little")
# 
# 
# medium.sample = sample(1:nrow(functional.community.matrix.relative.subsample),1000, replace = F)
# functional.community.matrix.relative.medium.subsample = functional.community.matrix.relative.subsample[medium.sample,]
# 
# system.time( # 165s
#   functional.community.matrix.relative.medium.subsample.SQchord <- distance(functional.community.matrix.relative.medium.subsample[,columns.with.functional.groups], method = "SQchord", weights = NULL, R = NULL, dist = T, fast = T))
# save(functional.community.matrix.relative.medium.subsample.SQchord,file = "./3_Output/IAa4_Calculate_Cell_Distances/IAa4_dist_matrix_abundance_rel_distance_SQchord_matrix_medium")
# 
# 
large.sample = sample(1:nrow(functional.community.matrix.relative.subsample),10000, replace = F)
functional.community.matrix.relative.large.subsample = functional.community.matrix.relative.subsample[large.sample,]

system.time( # 165s
  functional.community.matrix.relative.large.subsample.SQchord <- distance(functional.community.matrix.relative.large.subsample[,columns.with.functional.groups], method = "SQchord", weights = NULL, R = NULL, dist = T, fast = T))
save(functional.community.matrix.relative.large.subsample.SQchord,file = "./3_Output/IAa4_Calculate_Cell_Distances/IAa4_dist_matrix_abundance_rel_distance_SQchord_matrix_large")


## – End Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
end_time <- Sys.time()
end_time - start_time


