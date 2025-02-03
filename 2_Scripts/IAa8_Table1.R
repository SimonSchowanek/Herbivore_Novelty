######################################################################!
#                                                                    #
#      IAa10: MAKE  TABLE                                           #  
#                                                                    #
######################################################################!
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # METADATA
    # This script creates a table with summary stats about the functional types.
    #
    # Output: a table 
    # time:  less than one minute
    # Status: last update (2025/02/03) 
    # Note: 
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
# – Setting working directory ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
getwd()


# Macropus rufus has been renamed as Osphranter rufus 

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# – Setting up libraries ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# GIS
library("rworldmap") # visualisation data
library(raster) # GIS
library(sp) # GIS
library(rgeos) # GIS
library(rgdal) # GIS

# Misc
library(plyr) # to use JOIN function

# For the pretty table
library(formattable)
library(htmltools) # to save the pretty table
library(webshot) # to save the pretty table



# – Load Data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
functional.community.matrix.relative = read.csv("./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_relative_abundances.csv")
functional.community.matrix.relative = functional.community.matrix.relative[,!(colnames(functional.community.matrix.relative) %in% c("X","CellNr"))] # remove two redundant variables


functional.community.matrix.pres.abs = read.csv("./3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_presence_absence.csv")
functional.community.matrix.pres.abs = functional.community.matrix.pres.abs[,!(colnames(functional.community.matrix.pres.abs) %in% c("X","CellNr"))] # remove two redundant variables

herbivores.4 = read.csv("./3_Output/IAa2_Manual_Functional_Groups/IAa2_herbivores.csv")

raster.example = raster("./1_Input/PHYLACINE_V1_2_1/Ranges/Current/Abditomys_latidens.tif")


## – Directories  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
if(!dir.exists(paste0("./3_Output/IAa8_Funtional_Types_table"))){
  dir.create(paste0("./3_Output/IAa8_Funtional_Types_table"))}


# – Correct Biogeo Realm ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
functional.community.matrix.relative = functional.community.matrix.relative[which(functional.community.matrix.relative$Herbivores.Present == 1),] # Select cells with herbivores
functional.community.matrix.relative$Biogeo.realm = as.factor(as.character(functional.community.matrix.relative$Biogeo.realm))



# Select the functional group columns
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
columns.with.functional.groups = which(grepl("Functional.Group.", colnames(functional.community.matrix.relative)))

colnames(functional.community.matrix.relative) = gsub(pattern = "Functional.Group.",
                                                      replacement = "",
                                                      x = colnames(functional.community.matrix.relative))

colnames(functional.community.matrix.pres.abs) = gsub(pattern = "Functional.Group.",
                                                       replacement = "",
                                                       x = colnames(functional.community.matrix.pres.abs))


functional.groups.new.order = c(
  #Meso
  "semi.aquatic.mesoherbivore.grazer",
  "semi.aquatic.mesoherbivore.mixedfeeder",
  "semi.aquatic.mesoherbivore.browser",
  "terrestrial.mesoherbivore.grazer",
  "terrestrial.mesoherbivore.mixedfeeder",
  "terrestrial.mesoherbivore.browser", 
  "terrestrial.mesoherbivore.omnivore",
  "semi.arboreal.mesoherbivore.browser",
  "semi.arboreal.mesoherbivore.omnivore",
  "arboreal.mesoherbivore.browser", 
  "arboreal.mesoherbivore.omnivore",
  #Macro
  "semi.aquatic.macroherbivore.grazer",
  "semi.aquatic.macroherbivore.mixedfeeder",
  "semi.aquatic.macroherbivore.browser",
  "terrestrial.macroherbivore.grazer",
  "terrestrial.macroherbivore.mixedfeeder",
  "terrestrial.macroherbivore.browser",
  "terrestrial.macroherbivore.omnivore",
  "semi.arboreal.macroherbivore.browser",
  # Mega
  "semi.aquatic.megaherbivore.grazer",
  "semi.aquatic.megaherbivore.mixedfeeder",
  "terrestrial.megaherbivore.grazer",
  "terrestrial.megaherbivore.mixedfeeder",
  "terrestrial.megaherbivore.browser")




# Change column order
functional.community.matrix.relative <- functional.community.matrix.relative[, c(functional.groups.new.order, setdiff(names(functional.community.matrix.relative), functional.groups.new.order))]
names(functional.community.matrix.relative)
 
functional.community.matrix.pres.abs <- functional.community.matrix.pres.abs[, c(functional.groups.new.order, setdiff(names(functional.community.matrix.pres.abs), functional.groups.new.order))]
names(functional.community.matrix.pres.abs)




# Separate scenarios
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
df.current = functional.community.matrix.relative[which(functional.community.matrix.relative$Scenario == "Current"),]
df.pn = functional.community.matrix.relative[which(functional.community.matrix.relative$Scenario == "Present-Natural"),]





#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 2. FUNCTIONAL GROUP TABLE ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####


#Create the data frame
list.of.functional.groups = colnames(df.pn[,1:length(columns.with.functional.groups)])
fg.table = data.frame(matrix(NA,length(columns.with.functional.groups),5))
colnames(fg.table) = c("FG", "Example.Species", "Cells.PN", "Cells.Current")
fg.table$FG = list.of.functional.groups
fg.table$FG <- factor(fg.table$FG, levels = list.of.functional.groups)

# Calculate Species richeness of each group
functional.group.richness = table(herbivores.4$Functional.Group.Manual)
group.richness = data.frame(functional.group.richness)
colnames(group.richness) = c("FG", "SR")

# Merge tables
fg.table = join(fg.table, group.richness, by = "FG")

# Reorder Columns
fg.table = fg.table[,  c("FG","Cells.PN", "Cells.Current","SR", "Example.Species")]

# Fill in Range Sizes
df.pn.pres.abs = functional.community.matrix.pres.abs[which(functional.community.matrix.pres.abs$Scenario == "Present-Natural"),]
df.current.pres.abs = functional.community.matrix.pres.abs[which(functional.community.matrix.pres.abs$Scenario == "Current"),]

# Run for loop to calculate the % of cells occupied by the functional group
for (i in 1:nrow(fg.table)) {
  
  fg.table[i,"Cells.PN"] = round(100 * sum(df.pn.pres.abs[,which(colnames(df.pn.pres.abs) %in% fg.table[i,"FG"])]) / nrow(df.pn.pres.abs), 2)
  fg.table[i,"Cells.Current"] = round(100 * sum(df.current.pres.abs[,which(colnames(df.current.pres.abs) %in% fg.table[i,"FG"])]) / nrow(df.pn.pres.abs),2) # we standardise by the total number of cells containing herbivore in the present-natural
  
}


fg.table[which(fg.table$FG == "terrestrial.mesoherbivore.omnivore"),]


# Remove groups whose distribution == 0 (= semi.arboreal.mesoherbivore.mixedfeeder)
fg.table = fg.table[which(!((fg.table$Cells.PN == 0) & (fg.table$Cells.Current == 0))),]
  
# fg.table$Cells.PN - fg.table$Cells.Current

## make a species list of extant herbivores
extant.herbivores = herbivores.4[which((herbivores.4$Species.Status %in% c("Extant"))),]

for (i.FG in 1:nrow(fg.table)) {
  
  # list the species that belong to now extinct functional groups
  v.FG = fg.table[i.FG,"FG"]
  species.list = herbivores.4[which(as.character(herbivores.4$Functional.Group.Manual) == v.FG),"Binomial"] # select species belonging to functional type
  
  if(any(species.list %in% extant.herbivores$Binomial)){ 
    species.list.extant = species.list[which(species.list %in% extant.herbivores$Binomial)]
    random.species = sample(species.list.extant, 1)
  } else{
    random.species = sample(species.list, 1)
    random.species = paste0(as.character(random.species)," (†)")
  }
  
  fg.table[i.FG,"Example.Species"] = paste0(as.character(random.species))
  
}

# Clean the Names of the Table
fg.table$Example.Species = gsub(pattern = "_",
                                replacement = " ",
                                x = fg.table$Example.Species)

fg.table$FG = gsub(pattern = "semi.aquatic",
                   replacement = "semi-aquatic",
                   x = fg.table$FG)

fg.table$FG = gsub(pattern = "semi.arboreal",
                   replacement = "semi-arboreal",
                   x = fg.table$FG)

fg.table$FG = gsub(pattern = "\\.",
                   replacement = " ",
                   x = fg.table$FG)

fg.table$FG = gsub(pattern = "mixedfeeder",
                   replacement = "mixed feeder",
                   x = fg.table$FG)




## SHOW WHICH FUNCTIONAL TYPES ARE EXTINCT
m = which(fg.table$Cells.Current == 0)
fg.table[m, "FG"] = paste0(fg.table[m, "FG"], " (†)")




# Save the data

## Create Maps Directory if it doesn't exist
if(!dir.exists(paste0("./3_Output/IAa8_Funtional_Types_table/"))){
  dir.create(paste0("./3_Output/IAa8_Funtional_Types_table/"))}

write.csv(fg.table, "./3_Output/IAa8_Funtional_Types_table/IAa8_functional_group_table.csv", row.names = F)





#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 4. PLOT RANGE DIFFERENCES ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

## Plot Range Differences
fg.table.2 = fg.table
str(fg.table.2)
fg.table.2$FG = as.factor(fg.table.2$FG)
fg.table.2$range.diff = fg.table.2$Cells.PN - fg.table.2$Cells.Current
plot(fg.table.2$FG, fg.table.2$range.diff, horizontal = T, las = 2)

v.aquatic <- grep("aquatic", fg.table.2$FG, value = TRUE)      # Extract elements with match
v.terrestrial <- grep("terrestrial", fg.table.2$FG, value = TRUE)      # Extract elements with match
v.arboreal <- grep("arboreal", fg.table.2$FG, value = TRUE)      # Extract elements with match

## Assign Strate
fg.table.2[which(fg.table.2$FG %in% v.aquatic),"Stratum"] = "aquatic"
fg.table.2[which(fg.table.2$FG %in% v.terrestrial),"Stratum"] = "terrestrial"
fg.table.2[which(fg.table.2$FG %in% v.arboreal),"Stratum"] = "arboreal"
fg.table.2$Stratum = factor(fg.table.2$Stratum, levels = c("aquatic", "terrestrial", "arboreal")) #turn into factor

## Plot
png(paste0("./3_Output/IAa8_Funtional_Types_table/IAa8_range_diff_functional_types.png"))
boxplot(fg.table.2$range.diff ~ fg.table.2$Stratum, horizontal = T, ylab = NULL, xlab = "Difference in cells between present-natural and current ranges")
dev.off()







#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 3. TABLE USING FORMATABBLE##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

#Change table names
colnames(fg.table) = c("Functional Type", "Cells Present-Natural (%)", "Cells Current (%)", "Species Richness", "Example Species")

#Set a few color variables to make our table more visually appealing
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"


row.names(fg.table) <- NULL

## Table
table.formatabble = formattable(fg.table, 
                                
                                align =c("l","r","r","r", "l"), 
                                
                                list(`Indicator Name` = formatter(
                                  "span", style = ~ style(color = "grey",font.weight = "bold")),
                                  `Cells Present-Natural (%)`= color_bar(customGreen),
                                  `Cells Current (%)`= color_bar(customGreen),
                                  `Species Richness`= color_bar(customRed)
                                  
                                ))



## Create function to save the pretty table
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

#' Export a Formattable as PNG, PDF, or JPEG
#'
#' @param f A formattable.
#' @param file Export path with extension .png, .pdf, or .jpeg.
#' @param width Width specification of the html widget being exported.
#' @param height Height specification of the html widget being exported.
#' @param background Background color specification.
#' @param delay Time to wait before taking webshot, in seconds.
#'
#' @importFrom formattable as.htmlwidget
#' @importFrom htmltools html_print
#' @importFrom webshot webshot
#'
#' @export
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

## Export the Table
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
export_formattable(
  table.formatabble,
  file = "./3_Output/IAa8_Funtional_Types_table/IAa8_table_unedited.png",
  width = "100%",
  height = NULL,
  background = "white",
  delay = 0.2
)

export_formattable(
  table.formatabble,
  file = "./3_Output/IAa8_Funtional_Types_table/IAa8_table_unedited.pdf",
  width = "100%",
  height = NULL,
  background = "white",
  delay = 0.2
)


