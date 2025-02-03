######################################################################!
#                                                                    #
#      IAa4: Make functional groups                                  #
#                                                                    #  
#                                                                    #
######################################################################!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# METADATA
# This script selects and appropriate number of clusters and then classifies herbivores into said amount of clusters (functional groups)
#
# Output: a herbivores traits file that also includes which cluster they belong to
# time: needs less than 1 min to run.
# Status: Ready
# Note: I ran this script on 2025/02/03
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


## – Start Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
start_time <- Sys.time() # Save Time


# – Loading data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
herbivores.3 = read.csv("./3_Output/IAa1_Create_Community_Matrices/IAa1_herbivores_2.csv") # only added to add clustes to herbivore file


## – Change Data  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––

        # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  
        #   Calling the following species arboreal is almost certainly a mistake. I make this correction to save us some problems later.
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

matches = which(herbivores.3$Binomial == "Amblyrhiza_inundata")
herbivores.3[matches, "Habitat.Arboreal"] = 0
herbivores.3[matches, "Habitat.Categorical"] = "Terrestrial"



## – Directories  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
if(!dir.exists(paste0("./3_Output/IAa2_Manual_Functional_Groups/"))){
  dir.create(paste0("./3_Output/IAa2_Manual_Functional_Groups/"))}





#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 2. MAKE MANUAL GROUPS ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
colnames(herbivores.3)

# – Divide herbivores by stratum ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
semi.aquatic.herbivores = herbivores.3[which(herbivores.3$Habitat.Categorical == "Semi.Aquatic"),]
terrestrial.herbivores = herbivores.3[which(herbivores.3$Habitat.Categorical == "Terrestrial"),]
semi.arboreal.herbivores = herbivores.3[which(herbivores.3$Habitat.Categorical == "Semi.Arboreal"),]
arboreal.herbivores = herbivores.3[which(herbivores.3$Habitat.Categorical == "Arboreal"),]

# – Divide herbivores by body mass ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
semi.aquatic.megaherbivores = semi.aquatic.herbivores[which(semi.aquatic.herbivores$Mass.g >= 1*10^6),] # 8 species
semi.aquatic.macroherbivores = semi.aquatic.herbivores[which((semi.aquatic.herbivores$Mass.g < 1*10^6) & (semi.aquatic.herbivores$Mass.g >= 1*10^5)),] # 23 species
semi.aquatic.mesoherbivores = semi.aquatic.herbivores[which(semi.aquatic.herbivores$Mass.g < 1*10^5),] # 14 species

# nrow(semi.aquatic.megaherbivores)
# nrow(semi.aquatic.macroherbivores)
# nrow(semi.aquatic.mesoherbivores)

terrestrial.megaherbivores = terrestrial.herbivores[which(terrestrial.herbivores$Mass.g >= 1*10^6),] # 40 species
terrestrial.macroherbivores = terrestrial.herbivores[which((terrestrial.herbivores$Mass.g < 1*10^6) & (terrestrial.herbivores$Mass.g >= 1*10^5)),] # 122 species
terrestrial.mesoherbivores = terrestrial.herbivores[which(terrestrial.herbivores$Mass.g < 1*10^5),] # 236 species

# nrow(terrestrial.megaherbivores)
# nrow(terrestrial.macroherbivores)
# nrow(terrestrial.mesoherbivores)

#semi.arboreal.megaherbivores = semi.arboreal.herbivores[which(semi.arboreal.herbivores$Mass.g >= 1*10^6),] # 0 species
semi.arboreal.macroherbivores = semi.arboreal.herbivores[which((semi.arboreal.herbivores$Mass.g < 1*10^6) & (semi.arboreal.herbivores$Mass.g >= 1*10^5)),] # 3 species
semi.arboreal.mesoherbivores = semi.arboreal.herbivores[which(semi.arboreal.herbivores$Mass.g < 1*10^5),] # 33 species

# nrow(semi.arboreal.megaherbivores)
# nrow(semi.arboreal.macroherbivores)
# nrow(semi.arboreal.mesoherbivores)

#arboreal.megaherbivores = arboreal.herbivores[which(arboreal.herbivores$Mass.g >= 1*10^6),] # 0 species
#arboreal.macroherbivores = arboreal.herbivores[which((arboreal.herbivores$Mass.g < 1*10^6) & (arboreal.herbivores$Mass.g >= 1*10^5)),] # 0 species
arboreal.mesoherbivores = arboreal.herbivores[which(arboreal.herbivores$Mass.g < 1*10^5),] # 23 species

# – Divide herbivores by diet consumption ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Semiaquatic 
semi.aquatic.megaherbivore.grazers = semi.aquatic.megaherbivores[which(semi.aquatic.megaherbivores$Guild.w.Omnivory == "Grazer"),] # 4 species
semi.aquatic.megaherbivore.mixedfeeders = semi.aquatic.megaherbivores[which(semi.aquatic.megaherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 2 species
semi.aquatic.megaherbivore.browsers = semi.aquatic.megaherbivores[which(semi.aquatic.megaherbivores$Guild.w.Omnivory == "Browser"),] # 2 species
semi.aquatic.megaherbivore.omnivore = semi.aquatic.megaherbivores[which(semi.aquatic.megaherbivores$Guild.w.Omnivory == "Omnivore"),] # 0 species

    # nrow(semi.aquatic.megaherbivore.grazers)
    # nrow(semi.aquatic.megaherbivore.mixedfeeders)
    # nrow(semi.aquatic.megaherbivore.browsers)
    # nrow(semi.aquatic.megaherbivore.omnivore)

semi.aquatic.macroherbivore.grazers = semi.aquatic.macroherbivores[which(semi.aquatic.macroherbivores$Guild.w.Omnivory == "Grazer"),] # 3 species
semi.aquatic.macroherbivore.mixedfeeders = semi.aquatic.macroherbivores[which(semi.aquatic.macroherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 8 species
semi.aquatic.macroherbivore.browsers = semi.aquatic.macroherbivores[which(semi.aquatic.macroherbivores$Guild.w.Omnivory == "Browser"),] # 12 species
semi.aquatic.macroherbivore.omnivores = semi.aquatic.macroherbivores[which(semi.aquatic.macroherbivores$Guild.w.Omnivory == "Omnivore"),] # 0 species

    # nrow(semi.aquatic.macroherbivore.grazers)
    # nrow(semi.aquatic.macroherbivore.mixedfeeders)
    # nrow(semi.aquatic.macroherbivore.browsers)
    # nrow(semi.aquatic.macroherbivore.omnivores)

semi.aquatic.mesoherbivore.grazers = semi.aquatic.mesoherbivores[which(semi.aquatic.mesoherbivores$Guild.w.Omnivory == "Grazer"),] # 5 species
semi.aquatic.mesoherbivore.mixedfeeders = semi.aquatic.mesoherbivores[which(semi.aquatic.mesoherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 5 species
semi.aquatic.mesoherbivore.browsers = semi.aquatic.mesoherbivores[which(semi.aquatic.mesoherbivores$Guild.w.Omnivory == "Browser"),] # 4 species
semi.aquatic.mesoherbivore.omnivores = semi.aquatic.mesoherbivores[which(semi.aquatic.mesoherbivores$Guild.w.Omnivory == "Omnivore"),] # 0 species

    # nrow(semi.aquatic.mesoherbivore.grazers)
    # nrow(semi.aquatic.mesoherbivore.mixedfeeders)
    # nrow(semi.aquatic.mesoherbivore.browsers)
    # nrow(semi.aquatic.mesoherbivore.omnivores)

# Terrestrial
terrestrial.megaherbivore.grazers = terrestrial.megaherbivores[which(terrestrial.megaherbivores$Guild.w.Omnivory == "Grazer"),] # 8 species
terrestrial.megaherbivore.mixedfeeders = terrestrial.megaherbivores[which(terrestrial.megaherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 25 species
terrestrial.megaherbivore.browsers = terrestrial.megaherbivores[which(terrestrial.megaherbivores$Guild.w.Omnivory == "Browser"),] # 7 species
terrestrial.megaherbivore.omnivores = terrestrial.megaherbivores[which(terrestrial.megaherbivores$Guild.w.Omnivory == "Omnivore"),] # 0 species

    # nrow(terrestrial.megaherbivore.grazers)
    # nrow(terrestrial.megaherbivore.mixedfeeders)
    # nrow(terrestrial.megaherbivore.browsers)
    # nrow(terrestrial.megaherbivore.omnivores)

terrestrial.macroherbivore.grazers = terrestrial.macroherbivores[which(terrestrial.macroherbivores$Guild.w.Omnivory == "Grazer"),] # 28 species
terrestrial.macroherbivore.mixedfeeders = terrestrial.macroherbivores[which(terrestrial.macroherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 51 species
terrestrial.macroherbivore.browsers = terrestrial.macroherbivores[which(terrestrial.macroherbivores$Guild.w.Omnivory == "Browser"),] # 34 species
terrestrial.macroherbivore.omnivores = terrestrial.macroherbivores[which(terrestrial.macroherbivores$Guild.w.Omnivory == "Omnivore"),] # 9 species

    # nrow(terrestrial.macroherbivore.grazers)
    # nrow(terrestrial.macroherbivore.mixedfeeders)
    # nrow(terrestrial.macroherbivore.browsers)
    # nrow(terrestrial.macroherbivore.omnivores)

terrestrial.mesoherbivore.grazers = terrestrial.mesoherbivores[which(terrestrial.mesoherbivores$Guild.w.Omnivory == "Grazer"),] # 27 species
terrestrial.mesoherbivore.mixedfeeders = terrestrial.mesoherbivores[which(terrestrial.mesoherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 91 species
terrestrial.mesoherbivore.browsers = terrestrial.mesoherbivores[which(terrestrial.mesoherbivores$Guild.w.Omnivory == "Browser"),] # 92 species
terrestrial.mesoherbivore.omnivores = terrestrial.mesoherbivores[which(terrestrial.mesoherbivores$Guild.w.Omnivory == "Omnivore"),] # 26 species

    # nrow(terrestrial.mesoherbivore.grazers)
    # nrow(terrestrial.mesoherbivore.mixedfeeders)
    # nrow(terrestrial.mesoherbivore.browsers)
    # nrow(terrestrial.mesoherbivore.omnivores)

# Semi-Arboreal
semi.arboreal.macroherbivore.grazers = semi.arboreal.macroherbivores[which(semi.arboreal.macroherbivores$Guild.w.Omnivory == "Grazer"),] # 0 species
semi.arboreal.macroherbivore.mixedfeeders = semi.arboreal.macroherbivores[which(semi.arboreal.macroherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 0 species
semi.arboreal.macroherbivore.browsers = semi.arboreal.macroherbivores[which(semi.arboreal.macroherbivores$Guild.w.Omnivory == "Browser"),] # 3 species
semi.arboreal.macroherbivore.omnivore = semi.arboreal.macroherbivores[which(semi.arboreal.macroherbivores$Guild.w.Omnivory == "Omnivore"),] # 0 species

    # nrow(semi.arboreal.macroherbivore.grazers)
    # nrow(semi.arboreal.macroherbivore.mixedfeeders)
    # nrow(semi.arboreal.macroherbivore.browsers)
    # nrow(semi.arboreal.macroherbivore.omnivore)


semi.arboreal.mesoherbivore.grazers = semi.arboreal.mesoherbivores[which(semi.arboreal.mesoherbivores$Guild.w.Omnivory == "Grazer"),] # 0 species
semi.arboreal.mesoherbivore.mixedfeeders = semi.arboreal.mesoherbivores[which(semi.arboreal.mesoherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 1 species
semi.arboreal.mesoherbivore.browsers = semi.arboreal.mesoherbivores[which(semi.arboreal.mesoherbivores$Guild.w.Omnivory == "Browser"),] # 26 species
semi.arboreal.mesoherbivore.omnivores = semi.arboreal.mesoherbivores[which(semi.arboreal.mesoherbivores$Guild.w.Omnivory == "Omnivore"),] # 6 species

    # nrow(semi.arboreal.mesoherbivore.grazers)
    # nrow(semi.arboreal.mesoherbivore.mixedfeeders)
    # nrow(semi.arboreal.mesoherbivore.browsers)
    # nrow(semi.arboreal.mesoherbivore.omnivores)

# Arboreal
arboreal.mesoherbivore.grazers = arboreal.mesoherbivores[which(arboreal.mesoherbivores$Guild.w.Omnivory == "Grazer"),] # 0 species
arboreal.mesoherbivore.mixedfeeders = arboreal.mesoherbivores[which(arboreal.mesoherbivores$Guild.w.Omnivory == "Mixed Feeder"),] # 0 species
arboreal.mesoherbivore.browsers = arboreal.mesoherbivores[which(arboreal.mesoherbivores$Guild.w.Omnivory == "Browser"),] # 22 species
arboreal.mesoherbivore.omnivores = arboreal.mesoherbivores[which(arboreal.mesoherbivores$Guild.w.Omnivory == "Omnivore"),] # 1 species

    # nrow(arboreal.mesoherbivore.grazers)
    # nrow(arboreal.mesoherbivore.mixedfeeders)
    # nrow(arboreal.mesoherbivore.browsers)
    # nrow(arboreal.mesoherbivore.omnivores)


# – Assign manual classification ----
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

herbivores.3$Functional.Group.Manual = NA

herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.megaherbivore.grazers$Binomial),"Functional.Group.Manual"] = "semi.aquatic.megaherbivore.grazer"
herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.megaherbivore.mixedfeeders$Binomial),"Functional.Group.Manual"] = "semi.aquatic.megaherbivore.mixedfeeder"
herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.megaherbivore.browsers$Binomial),"Functional.Group.Manual"] = "semi.aquatic.macroherbivore.browser"

herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.macroherbivore.grazers$Binomial),"Functional.Group.Manual"] = "semi.aquatic.macroherbivore.grazer"
herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.macroherbivore.mixedfeeders$Binomial),"Functional.Group.Manual"] = "semi.aquatic.macroherbivore.mixedfeeder"
herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.macroherbivore.browsers$Binomial),"Functional.Group.Manual"] = "semi.aquatic.macroherbivore.browser"

herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.mesoherbivore.grazers$Binomial),"Functional.Group.Manual"] = "semi.aquatic.mesoherbivore.grazer"
herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.mesoherbivore.mixedfeeders$Binomial),"Functional.Group.Manual"] = "semi.aquatic.mesoherbivore.mixedfeeder"
herbivores.3[which(herbivores.3$Binomial %in% semi.aquatic.mesoherbivore.browsers$Binomial),"Functional.Group.Manual"] = "semi.aquatic.mesoherbivore.browser"

herbivores.3[which(herbivores.3$Binomial %in% terrestrial.megaherbivore.grazers$Binomial),"Functional.Group.Manual"] = "terrestrial.megaherbivore.grazer"
herbivores.3[which(herbivores.3$Binomial %in% terrestrial.megaherbivore.mixedfeeders$Binomial),"Functional.Group.Manual"] = "terrestrial.megaherbivore.mixedfeeder"
herbivores.3[which(herbivores.3$Binomial %in% terrestrial.megaherbivore.browsers$Binomial),"Functional.Group.Manual"] = "terrestrial.megaherbivore.browser"

herbivores.3[which(herbivores.3$Binomial %in% terrestrial.macroherbivore.grazers$Binomial),"Functional.Group.Manual"] = "terrestrial.macroherbivore.grazer"
herbivores.3[which(herbivores.3$Binomial %in% terrestrial.macroherbivore.mixedfeeders$Binomial),"Functional.Group.Manual"] = "terrestrial.macroherbivore.mixedfeeder"
herbivores.3[which(herbivores.3$Binomial %in% terrestrial.macroherbivore.browsers$Binomial),"Functional.Group.Manual"] = "terrestrial.macroherbivore.browser"
herbivores.3[which(herbivores.3$Binomial %in% terrestrial.macroherbivore.omnivores$Binomial),"Functional.Group.Manual"] = "terrestrial.macroherbivore.omnivore"

herbivores.3[which(herbivores.3$Binomial %in% terrestrial.mesoherbivore.grazers$Binomial),"Functional.Group.Manual"] = "terrestrial.mesoherbivore.grazer"
herbivores.3[which(herbivores.3$Binomial %in% terrestrial.mesoherbivore.mixedfeeders$Binomial),"Functional.Group.Manual"] = "terrestrial.mesoherbivore.mixedfeeder"
herbivores.3[which(herbivores.3$Binomial %in% terrestrial.mesoherbivore.browsers$Binomial),"Functional.Group.Manual"] = "terrestrial.mesoherbivore.browser"
herbivores.3[which(herbivores.3$Binomial %in% terrestrial.mesoherbivore.omnivores$Binomial),"Functional.Group.Manual"] = "terrestrial.mesoherbivore.omnivore"

herbivores.3[which(herbivores.3$Binomial %in% semi.arboreal.macroherbivore.browsers$Binomial),"Functional.Group.Manual"] = "semi.arboreal.macroherbivore.browser"

herbivores.3[which(herbivores.3$Binomial %in% semi.arboreal.mesoherbivore.mixedfeeders$Binomial),"Functional.Group.Manual"] = "semi.arboreal.mesoherbivore.mixedfeeder"
herbivores.3[which(herbivores.3$Binomial %in% semi.arboreal.mesoherbivore.browsers$Binomial),"Functional.Group.Manual"] = "semi.arboreal.mesoherbivore.browser"
herbivores.3[which(herbivores.3$Binomial %in% semi.arboreal.mesoherbivore.omnivores$Binomial),"Functional.Group.Manual"] = "semi.arboreal.mesoherbivore.omnivore"

herbivores.3[which(herbivores.3$Binomial %in% arboreal.mesoherbivore.browsers$Binomial),"Functional.Group.Manual"] = "arboreal.mesoherbivore.browser"
herbivores.3[which(herbivores.3$Binomial %in% arboreal.mesoherbivore.omnivores$Binomial),"Functional.Group.Manual"] = "arboreal.mesoherbivore.omnivore"


any(is.na(herbivores.3$Functional.Group.Manual)) # No NA's left
length(unique(herbivores.3$Functional.Group.Manual))

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 3. SAVE GROUPS ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

write.csv(herbivores.3, "./3_Output/IAa2_Manual_Functional_Groups/IAa2_herbivores.csv", row.names=FALSE)

## – End Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
end_time <- Sys.time()
end_time - start_time
