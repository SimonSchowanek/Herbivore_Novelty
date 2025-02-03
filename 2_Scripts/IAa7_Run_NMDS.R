######################################################################!
#                                                                    #
#      IAa7: Plot communities with similar structure                 #  
#                                                                    #
######################################################################!


      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # METADATA
      # This script performs a NMDS and saves the resulting output (make sure the outrepo is empty)
      #
      # Output: a csv file with scores and a NMDS object (and a plot)
      # time:  takes a while to run (10-15 min)
      # Status: Ready (2025/03/03)
      # Notes: 
      #
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
#setwd("/Users/au572919/Research/I_MegaPast2Future/IA_Herbivore_Functional_Groups/IA_Herbivore_Functional_Groups/")


## – Start Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
start_time <- Sys.time() # Save Time



# – Setting up libraries ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr) # for mutiple panels
library(ggrepel) # for cleaner labels on the graph



# – Load Data ####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
load("./3_Output/IAa4_Calculate_Cell_Distances/IAa4_dist_matrix_abundance_rel_distance_SQchord")

functional.community.matrix = read.csv("3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_relative_abundances_weighed.csv")
colnames(functional.community.matrix)
functional.community.matrix$ID = rownames(functional.community.matrix) # add unique identifier
matches = which(!(grepl("Functional.Group.", colnames(functional.community.matrix)))) # select columns that are not functional groups (i.e. metadata)
metadata.functional.community.matrix = functional.community.matrix[,matches]

functional.community.matrix.absolute = read.csv("3_Output/IAa3_Create_Functional_Community_Matrix/IAa3_functional_community_matrix_abundances.csv")


## – Set parameters ####
#–––––––––––––––––––––––––––––––––––––––––––––––––

## Visual Parameters
theme_set(theme_bw()) # set ggplot layout
p.nr.functional.groups = length(which((grepl("Functional.Group.", colnames(functional.community.matrix))))) # select columns that are not functional groups (i.e. metadata)



## – Correct data ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
functional.community.matrix[which(functional.community.matrix$Biogeo.realm == "Africotropical"),"Biogeo.realm"] = "Afrotropical"



## – Directories  ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
if(!dir.exists(paste0("./3_Output/IAa7_Run_NMDS/"))){
  dir.create(paste0("./3_Output/IAa7_Run_NMDS/"))}




#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 2. CALCULATE FUNCTIONAL RICHNESS ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

## Select columns with functional groups
matches = (which((grepl("Functional.Group.", colnames(functional.community.matrix.absolute))))) # select columns that are not functional groups (i.e. metadata)

## Change group richness to presence/absence data for functional groups
df.temp = functional.community.matrix.absolute[,matches]
df.temp[df.temp > 0] = 1 # Turn all values higher than O into ones

## Calculate Richness of Functional Types
functional.community.matrix.absolute$functional.type.richness = apply(df.temp,1,sum)
str(functional.community.matrix.absolute)
rm(df.temp)

# Make Violin Plot
set.seed(123) # to save the jitter plot from moving



(p.FRic <- ggplot(functional.community.matrix.absolute[functional.community.matrix.absolute$Herbivores.Present == 1,], aes(x=Scenario, y=functional.type.richness, fill=Scenario)) + # fill=name allow to automatically dedicate a color for each group
    geom_jitter(alpha = 0.1, width = .2) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          legend.position = "none") +
    xlab("Richness of Herbivore Functional Types") +
    coord_flip() +
    geom_boxplot(alpha = 0.8, fill = "red", col = "#8b0000", width = 0.2))



(p.FRic <- ggplot(data = functional.community.matrix.absolute[functional.community.matrix.absolute$Herbivores.Present == 1,]) + # fill=name allow to automatically dedicate a color for each group
    geom_histogram(data = functional.community.matrix.absolute[functional.community.matrix.absolute$Herbivores.Present == 1,],
                   aes(x=functional.type.richness), binwidth = 1, colour = "black") + 
    theme(axis.title.y = element_blank(),
          legend.position = "none",
          strip.text = element_text(size=18, face = "bold"),
          axis.title.x = element_text(size=15)) +
    xlab("Richness of Herbivore Functional Types") +
    facet_wrap(vars(Scenario), ncol = 1))


#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 3. RUN THE NMDS ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

## Select cells where herbivores are present
functional.community.matrix = functional.community.matrix[which(functional.community.matrix$Herbivores.Present == 1),]

## divide the df into the two scenarios and select 500 samples in each scenario
functional.community.matrix.Current = functional.community.matrix[which(functional.community.matrix$Scenario == "Current"),]
functional.community.matrix.PN = functional.community.matrix[which(functional.community.matrix$Scenario == "Present-Natural"),]
samples.current = sample(1:nrow(functional.community.matrix.Current), 1000) # sample
samples.PN = sample(1:nrow(functional.community.matrix.PN), 1000) # sample
functional.community.matrix.subset.current = functional.community.matrix.Current[samples.current, 1:p.nr.functional.groups] # Select 500 samples
functional.community.matrix.subset.PN = functional.community.matrix.PN[samples.PN,1:p.nr.functional.groups] # select 500 samples
functional.community.matrix.subset = rbind(functional.community.matrix.subset.current,functional.community.matrix.subset.PN) # bind the selected samples back together

####
system.time(example_NMDS <- metaMDS(functional.community.matrix.subset, # Our community-by-species matrix
                                    k=3,
                                    trymax= 100, 
                                    autotransform = T)) # The number of reduced dimensions







#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####
##### 4. PLOT THE NMDS ##### 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––####

# exploratory approach

stressplot(example_NMDS)
# Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions


#plot(example_NMDS)
# The plot shows us both the communities (“sites”, open circles) and functional groups (red crosses), but we don’t know which circle corresponds to which site, and which species corresponds to which cross.
# We can use the function ordiplot and orditorp to add text to the plot in place of points to make some sense of this rather non-intuitive mess.


## Collect Point Data
df.points = data.frame(example_NMDS$points)
tail(df.points)
df.points$label = as.factor(rownames(df.points))
df.points = merge(df.points, metadata.functional.community.matrix, by.x = "label", by.y = "ID")

## Collect Point Data
df.species = data.frame(example_NMDS$species)
tail(df.species)
df.species$functional.group = as.factor(rownames(df.species))
df.species$functional.group = gsub(x = df.species$functional.group, pattern = "Functional.Group.", replacement = "") # Remove redundant text
df.species$functional.group.2 = NA


## Dive into size classes
matches = which((grepl("mega", df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Size.class"] = "megaherbivore"
matches = which((grepl("macro",  df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Size.class"] = "macroherbivore"
matches = which((grepl("meso",  df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Size.class"] = "mesoherbivore"

## Divide into diet classes
matches = which((grepl("grazer", df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Diet.class"] = "grazer"
matches = which((grepl("mixed",  df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Diet.class"] = "mixed feeder"
matches = which((grepl("browser",  df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Diet.class"] = "browser"
matches = which((grepl("omnivore",  df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Diet.class"] = "omnivore"

## Divide into substrate
matches = which((grepl("semi.aquatic", df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Habitat.class"] = "semi.aquatic"
matches = which((grepl("terrestrial",  df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Habitat.class"] = "terrestrial"
matches = which((grepl("semi.arboreal",  df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Habitat.class"] = "semi.arboreal"
matches = which((grepl("arboreal",  df.species$functional.group))) # select columns that are not functional groups (i.e. metadata)
df.species[matches, "Habitat.class"] = "arboreal"


v.herbivore.sequence = c("semi.aquatic.mesoherbivore.grazer",
                         "semi.aquatic.mesoherbivore.mixedfeeder",
                         "semi.aquatic.mesoherbivore.browser",
                         "terrestrial.mesoherbivore.grazer", 
                         "terrestrial.mesoherbivore.mixedfeeder",
                         "terrestrial.mesoherbivore.browser",
                         "terrestrial.mesoherbivore.omnivore",
                         #"semi.arboreal.mesoherbivore.mixedfeeder",
                         "semi.arboreal.mesoherbivore.browser",
                         "semi.arboreal.mesoherbivore.omnivore",
                         "arboreal.mesoherbivore.browser",
                         "arboreal.mesoherbivore.omnivore",
                         "semi.aquatic.macroherbivore.grazer",
                         "semi.aquatic.macroherbivore.mixedfeeder",
                         "semi.aquatic.macroherbivore.browser",
                         "terrestrial.macroherbivore.grazer",
                         "terrestrial.macroherbivore.mixedfeeder",
                         "terrestrial.macroherbivore.browser",
                         "terrestrial.macroherbivore.omnivore",
                         "semi.arboreal.macroherbivore.browser",
                         "semi.aquatic.megaherbivore.grazer",
                         "semi.aquatic.megaherbivore.mixedfeeder",
                         "terrestrial.megaherbivore.grazer",
                         "terrestrial.megaherbivore.mixedfeeder",
                         "terrestrial.megaherbivore.browser")


df.species$functional.group = factor(df.species$functional.group, levels = v.herbivore.sequence) # change factor sequence
df.species  = df.species[order(df.species$functional.group),] # change factor sequence
df.species$ID = 1:nrow(df.species)

## Select only the most influential species (n = 10)
df.species$vector.length.3D = sqrt(df.species$MDS1^2 + df.species$MDS2^2 + df.species$MDS3^2)
df.species = df.species[order(df.species$vector.length.3D, decreasing = T),]
write.csv(df.species, file = "./3_Output/IAa7_Run_NMDS/IAa7_NMDS_Loadings.csv")

#df.species = df.species[1:5,]

    # ~~~ EXPLANATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  
    #   We create the same plot twice, but we set alpha to 0 for the points and the vectors respectively, to make sure the dimensions in both plots are the same.
    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## – Plot Parameters data ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
title.position.x = -0.2
title.position.y = 2
gridd.cross.width = 0.5
loadings.alpha = 0.7
loadings.width = 1.2
label.size = 4



## – Calculate Hulls Parameters data ####
#––––––––––––––––––––––––––––––––––––––––––––––––– 
hull.PN <- subset(df.points, Scenario=="Present-Natural") %>% 
  slice(chull(MDS1, MDS2))

hull.Current <- subset(df.points, Scenario=="Current") %>% 
  slice(chull(MDS1, MDS2))


# Function to compute convex hull per Biogeographic Realm
hull_df.PN <- df.points %>%
  group_by(Biogeo.realm) %>%
  filter(Scenario == "Present-Natural") %>%
  slice(chull(MDS1, MDS2))  # Compute convex hull

# Function to compute convex hull per Biogeographic Realm
hull_df.Current <- df.points %>%
  group_by(Biogeo.realm) %>%
  filter(Scenario == "Current") %>%
  slice(chull(MDS1, MDS2))  # Compute convex hull






## – SUMMARY POINTS ####
#–––––––––––––––––––––––––––––––––––––––––––––––––

# PN
points.PN = subset(df.points, Scenario=="Present-Natural")

df.summary.points.PN = data.frame(
  MDS1 = median(points.PN$MDS1),
  MDS2 = median(points.PN$MDS2))

df.summary.points.PN$MDS1_Q2.5 = quantile(points.PN$MDS1, c(0.025))
df.summary.points.PN$MDS1_Q97.5 = quantile(points.PN$MDS1, c(0.975))

df.summary.points.PN$MDS2_Q2.5 = quantile(points.PN$MDS2, c(0.025))
df.summary.points.PN$MDS2_Q97.5 = quantile(points.PN$MDS2, c(0.975))




# CURRENT
points.current = subset(df.points, Scenario=="Current")

df.summary.points.Current = data.frame(
  MDS1 = median(points.current$MDS1),
  MDS2 = median(points.current$MDS2))

df.summary.points.Current$MDS1_Q2.5 = quantile(points.current$MDS1, c(0.025))
df.summary.points.Current$MDS1_Q97.5 = quantile(points.current$MDS1, c(0.975))

df.summary.points.Current$MDS2_Q2.5 = quantile(points.current$MDS2, c(0.025))
df.summary.points.Current$MDS2_Q97.5 = quantile(points.current$MDS2, c(0.975))






## – Plot Panels ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
legend.text.size = 6

(p1 = ggplot() +
    coord_fixed() +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = gridd.cross.width, colour  = "black") +
    geom_vline(aes(xintercept = 0), linetype = "dashed", size = gridd.cross.width, colour  = "black") +
    
    # geom_polygon(data = hull.PN, alpha = 0.2, # a hull around all points
    #                 aes(x = MDS1, y = MDS2), fill = "grey85", colour = "black",
    #               linetype = "dashed") +
    
    geom_point(data = subset(df.points, Scenario=="Current"),
               aes(x= MDS1, y = MDS2, fill = Biogeo.realm,
                   col = Biogeo.realm),
               alpha = 0, 
               shape = 6,
               size = 3) + 
    
    geom_point(data = subset(df.points, Scenario=="Present-Natural"),
               aes(x= MDS1, y = MDS2, fill = Biogeo.realm,
                   col = Biogeo.realm),
               alpha = 1, 
               shape = 1,
               size = 3) +
    
    geom_polygon(data = hull_df.PN, aes(x = MDS1, y = MDS2, fill = Biogeo.realm, colour = Biogeo.realm, group = Biogeo.realm),
                 alpha = 0, linetype = "dashed", size = 0.8) +
    
    geom_point(data = df.summary.points.PN,
               aes(x= MDS1,
                   y = MDS2),
               fill = "black",
               colour = "black",
               alpha = 1, 
               shape = 16,
               size = 6) +
    
    geom_segment(data = df.summary.points.PN,
                 mapping = aes(x= MDS1_Q2.5, y = MDS2, xend = MDS1_Q97.5, yend = MDS2),
                 colour = "black",
                 size = .8,
                 alpha = 1) +
    
    geom_segment(data = df.summary.points.PN,
                 mapping = aes(y= MDS2_Q2.5, x = MDS1, yend = MDS2_Q97.5, xend = MDS1),
                 colour = "black",
                 size = .8,
                 alpha = 1) +
    
    # geom_segment(data = df.species,
    #              mapping = aes(x= 0, y =0, xend= MDS1, yend = MDS2),
    #              arrow = arrow (length = unit(0.015, "npc"),
    #                             type = "closed"),
    #              colour = "black",
    #              size = 0,
    #              alpha = 0) +
    
    # geom_text(data = df.species,
    #           mapping = aes(label= ID, x = MDS1*1.1, MDS2 *1.1),
    #           colour = "black",
    #           alpha = 0) +
    
    annotate("text",
             -Inf, Inf,
             label = "Present-Natural",
             hjust = title.position.x+0.05,
             vjust = title.position.y,
             size = legend.text.size,
             fontface = "bold"))

# annotate("text",
#          -Inf, Inf,
#          label = paste0("Stress: ", round(example_NMDS$stress,2)),
#          hjust = title.position.x,
#          vjust = title.position.y + 2.5,
#          size = 5,
#          fontface = "italic"))



## Save legend  
p1_legend <- get_legend(p1)


(p2 = ggplot() +
    coord_fixed() +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = gridd.cross.width, colour  = "black") +
    geom_vline(aes(xintercept = 0), linetype = "dashed", size = gridd.cross.width, colour  = "black") +
    
    # geom_polygon(data = hull.Current, alpha = 0.2, # a hull around all points
    #              aes(x = MDS1, y = MDS2),
    #              fill = "grey85", colour = "black",
    #              linetype = "dashed") +
    #   
    geom_point(data = subset(df.points, Scenario=="Current"),
               aes(x= MDS1, y = MDS2, fill = Biogeo.realm,
                   col = Biogeo.realm),
               alpha = 0.8, 
               shape = 1,
               size = 3) + 
    
    geom_point(data = subset(df.points, Scenario=="Present-Natural"),
               aes(x= MDS1, y = MDS2, fill = Biogeo.realm,
                   col = Biogeo.realm),
               alpha = 0, 
               shape = 16,
               size = 3) +
    
    geom_polygon(data = hull_df.Current, aes(x = MDS1, y = MDS2, fill = Biogeo.realm, colour = Biogeo.realm, group = Biogeo.realm),
                 alpha = 0, linetype = "dashed", size = 0.8) +
    
    geom_point(data = df.summary.points.Current,
               aes(x= MDS1,
                   y = MDS2),
               fill = "black",
               colour = "black",
               alpha = 1, 
               shape = 16,
               size = 6) +
    
    geom_segment(data = df.summary.points.Current,
                 mapping = aes(x= MDS1_Q2.5, y = MDS2, xend = MDS1_Q97.5, yend = MDS2),
                 colour = "black",
                 size = 0.8,
                 alpha = 1) +
    
    geom_segment(data = df.summary.points.Current,
                 mapping = aes(y= MDS2_Q2.5, x = MDS1, yend = MDS2_Q97.5, xend = MDS1),
                 colour = "black",
                 size = 0.8,
                 alpha = 1) +
    
    theme(legend.position = c(0.2, 0.25),
          legend.background = element_rect(fill="white",
                                           colour ="black",
                                           size=0.5, linetype="solid")) +
    #   
    # geom_segment(data = df.species,
    #              mapping = aes(x= 0, y =0, xend= MDS1, yend = MDS2),
    #              arrow = arrow (length = unit(0.015, "npc"),
    #                             type = "closed"),
    #              colour = "black",
    #              size = 0,
    #              alpha = 0) + 
    #   
    # geom_text(data = df.species,
    #           mapping = aes(label= ID, x = MDS1*1.1, MDS2 *1.1),
    #           colour = "black",
    #           alpha = 0) +
    
    annotate("text",
             -Inf, Inf,
             label = "Current",
             hjust = title.position.x,
             vjust = title.position.y,
             size = legend.text.size,
             fontface = "bold"))







#- LOADINGS





## Body Sizes
l.temp = split(df.species, f = df.species$Size.class)
l.temp = lapply(l.temp, function(x){ apply(x[,c("MDS1", "MDS2", "MDS3")], 2, mean) })
df.loadings1 = do.call(rbind, l.temp)

## Diet
l.temp = split(df.species, f = df.species$Diet.class)
l.temp = lapply(l.temp, function(x){ apply(x[,c("MDS1", "MDS2", "MDS3")], 2, mean) })
df.loadings2 = do.call(rbind, l.temp)


## Terrestriality
l.temp = split(df.species, f = df.species$Habitat.class)
l.temp = lapply(l.temp, function(x){ apply(x[,c("MDS1", "MDS2", "MDS3")], 2, mean) })
df.loadings3 = do.call(rbind, l.temp)


df.loadings = rbind(df.loadings1, df.loadings2, df.loadings3)
df.loadings = data.frame(df.loadings)
df.loadings$loading = rownames(df.loadings)



df.species.subset = df.species


(p3 = ggplot() +
    coord_fixed() +
    
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = gridd.cross.width, colour  = "black") +
    
    geom_vline(aes(xintercept = 0), linetype = "dashed", size = gridd.cross.width, colour  = "black") +
    
    geom_point(data = subset(df.points, Scenario=="Current"),
               aes(x= MDS1, y = MDS2, fill = Biogeo.realm,
                   col = Biogeo.realm),
               alpha = 0, 
               shape = 6,
               size = 2.5) + 
    
    geom_point(data = subset(df.points, Scenario=="Present-Natural"),
               aes(x= MDS1, y = MDS2, fill = Biogeo.realm,
                   col = Biogeo.realm),
               alpha = 0, 
               shape = 16,
               size = 2.5) +
    
    geom_segment(data = df.loadings,
                 mapping = aes(x= 0, y =0, xend= MDS1, yend = MDS2),
                 arrow = arrow(length = unit(0.015, "npc"),
                               type = "closed"),
                 colour = "black",
                 size = loadings.width,
                 alpha = loadings.alpha,
                 linetype = "solid") +
    
    geom_label_repel(data = df.loadings,
                     mapping = aes(label= loading, x = MDS1, MDS2),
                     colour = "black",
                     alpha = 0.8,
                     nudge_x = 0,
                     nudge_y = 0,
                     size = label.size, 
                     box.padding = 0.4,
                     direction = "both",
                     max.overlaps = p.nr.functional.groups) +
    
    theme(legend.position="none") +
    
    annotate("text",
             -Inf, Inf,
             label = "Loadings",
             hjust = title.position.x,
             vjust = title.position.y,
             size = legend.text.size,
             fontface = "bold"))





# Plot Legend
p4 = ggplot() + geom_point(data = subset(df.points, Scenario=="Current"),
                           aes(x= MDS1, y = MDS2, fill = Biogeo.realm,
                               col = Biogeo.realm),
                           alpha = 0,
                           shape = 6,
                           size = 1.5) +
  geom_point(data = subset(df.points, Scenario=="Present-Natural"),
             aes(x= MDS1, y = MDS2, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0,
             shape = 16,
             size = 2.5) +
  # annotate(geom = "table", x = min(df.points$MDS1), y = max(df.points$MDS2), label = list(df.table),
  #                    vjust = 1, hjust = 0, size = 3) +
  theme_void() +
  theme(legend.position="none")


# get_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}

legend = ggarrange(p4,p1_legend)

(plot.final.1 = ggarrange(p2,
                          p1 + theme(legend.position="none"),
                          p3 + theme(legend.position="none"),
                          #legend,
                          ncol = 1,
                          nrow = 3,
                          labels = c("A", "B", "C")))










## – Plot Panels ####
#–––––––––––––––––––––––––––––––––––––––––––––––––






(plot.final.4 = ggarrange(plot.final.1,
                          p.FRic,
                          ncol = 2,
                          nrow = 1,
                          heights = c(3, 1),
                          labels = c(NA,"D")))



png(width=1000,height=1000, file = "3_Output/IAa7_Run_NMDS/IAa7_NMDS_Output_dim1_2_FRIC.png")
plot.final.4
dev.off()



pdf(width=15,height=15, file = "3_Output/IAa7_Run_NMDS/IAa7_NMDS_Output_dim1_2_FRIC.pdf")
plot.final.4
dev.off()






p3 = ggplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = gridd.cross.width, colour  = "black") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = gridd.cross.width, colour  = "black") +
  geom_point(data = subset(df.points, Scenario=="Current"),
             aes(x= MDS1, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 6,
             size = 2.5) + 
  geom_point(data = subset(df.points, Scenario=="Present-Natural"),
             aes(x= MDS1, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 16,
             size = 2.5) +
  geom_segment(data = subset(df.species, Size.class=="mesoherbivore"),
               mapping = aes(x= 0, y =0, xend= MDS1, yend = MDS3),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "darkorange",
               size = loadings.width,
               alpha = loadings.alpha,
               linetype = "solid") +
  geom_segment(data = subset(df.species, Size.class=="macroherbivore"),
               mapping = aes(x= 0, y =0, xend= MDS1, yend = MDS3),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "darkred",
               size = loadings.width,
               alpha = loadings.alpha,
               linetype = "solid") +
  geom_segment(data = subset(df.species, Size.class=="megaherbivore"),
               mapping = aes(x= 0, y =0, xend= MDS1, yend = MDS3),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "black",
               size = loadings.width,
               alpha = loadings.alpha,
               linetype = "solid") + 
  geom_label_repel(data = df.species,
                   mapping = aes(label= functional.group, x = MDS1, MDS3),
                   colour = "black",
                   alpha = 0.8,
                   nudge_x = 0,
                   nudge_y = 0,
                   size = label.size, 
                   box.padding = 0.5,
                   direction = "both",
                   max.overlaps = p.nr.functional.groups) +
  
  theme(legend.position="none") +
  
  annotate("text",
           -Inf, Inf,
           label = "Loadings",
           hjust = title.position.x,
           vjust = title.position.y,
           size = 5,
           fontface = "bold")


df.table = df.species[,c("ID", "functional.group")]

## Plot Legend
p4 = ggplot() + 
  
  geom_point(data = subset(df.points, Scenario=="Current"),
             aes(x= MDS1, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 6,
             size = 1.5) +
  
  geom_point(data = subset(df.points, Scenario=="Present-Natural"),
             aes(x= MDS1, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 16,
             size = 2.5) +
  
  # geom_table((geom = "text", x = min(df.points$MDS1), y = max(df.points$MDS3), label = df.table,
  #          vjust = 1, hjust = 0, size = 3) +
  
  theme_void() +
  theme(legend.position="none")


# get_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}

legend = ggarrange(p4,p1_legend)

(plot.final.2 = ggarrange(p1 + theme(legend.position="none"),
                          p2 + theme(legend.position="none"),
                          p3 + theme(legend.position="none"),
                          legend,
                          ncol = 2,
                          nrow = 2,
                          labels = c("A", "B", "C")))



png(width=1000,height=1000, file = "3_Output/IAa7_Run_NMDS/IAa7_NMDS_Output_dim1_3.png")
plot.final.2
dev.off()


pdf(width=10,height=10, file = "3_Output/IAa7_Run_NMDS/IAa7_NMDS_Output_dim1_3.pdf")
plot.final.2
dev.off()








## – Plot Panels ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
p1 = ggplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = gridd.cross.width, colour  = "black") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = gridd.cross.width, colour  = "black") +
  geom_point(data = subset(df.points, Scenario=="Current"),
             aes(x= MDS2, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 6,
             size = 2.5) + 
  geom_point(data = subset(df.points, Scenario=="Present-Natural"),
             aes(x= MDS2, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 1, 
             shape = 16,
             size = 2.5) +
  geom_segment(data = df.species,
               mapping = aes(x= 0, y =0, xend= MDS2, yend = MDS3),
               arrow = arrow (length = unit(0.015, "npc"),
                              type = "closed"),
               colour = "black",
               size = 0,
               alpha = 0) + 
  geom_text(data = df.species,
            mapping = aes(label= ID, x = MDS2*1.1, MDS3 *1.1),
            colour = "black",
            alpha = 0) +
  annotate("text",
           -Inf, Inf,
           label = "Present-Natural",
           hjust = title.position.x+0.05,
           vjust = title.position.y,
           size = 5,
           fontface = "bold") +
  annotate("text",
           -Inf, Inf,
           label = paste0("Stress: ", round(example_NMDS$stress,2)),
           hjust = title.position.x,
           vjust = title.position.y + 2.5,
           size = 5,
           fontface = "italic")

## Save legend  
p1_legend <- get_legend(p1)


p2 = ggplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = gridd.cross.width, colour  = "black") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = gridd.cross.width, colour  = "black") +
  geom_point(data = subset(df.points, Scenario=="Current"),
             aes(x= MDS2, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 1, 
             shape = 17,
             size = 2.5) + 
  geom_point(data = subset(df.points, Scenario=="Present-Natural"),
             aes(x= MDS2, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 16,
             size = 2.5) +
  geom_segment(data = df.species,
               mapping = aes(x= 0, y =0, xend= MDS2, yend = MDS3),
               arrow = arrow (length = unit(0.015, "npc"),
                              type = "closed"),
               colour = "black",
               size = 0,
               alpha = 0) + 
  geom_text(data = df.species,
            mapping = aes(label= ID, x = MDS2*1.1, MDS3 *1.1),
            colour = "black",
            alpha = 0) +
  theme(legend.title=element_blank()) +
  annotate("text",
           -Inf, Inf,
           label = "Current",
           hjust = title.position.x,
           vjust = title.position.y,
           size = 5,
           fontface = "bold")

p3 = ggplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = gridd.cross.width, colour  = "black") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = gridd.cross.width, colour  = "black") +
  geom_point(data = subset(df.points, Scenario=="Current"),
             aes(x= MDS2, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 6,
             size = 2.5) + 
  geom_point(data = subset(df.points, Scenario=="Present-Natural"),
             aes(x= MDS2, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 16,
             size = 2.5) +
  geom_segment(data = subset(df.species, Size.class=="mesoherbivore"),
               mapping = aes(x= 0, y =0, xend= MDS2, yend = MDS3),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "darkorange",
               size = loadings.width,
               alpha = loadings.alpha,
               linetype = "solid") +
  geom_segment(data = subset(df.species, Size.class=="macroherbivore"),
               mapping = aes(x= 0, y =0, xend= MDS2, yend = MDS3),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "darkred",
               size = loadings.width,
               alpha = loadings.alpha,
               linetype = "solid") +
  geom_segment(data = subset(df.species, Size.class=="megaherbivore"),
               mapping = aes(x= 0, y =0, xend= MDS2, yend = MDS3),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "black",
               size = loadings.width,
               alpha = loadings.alpha,
               linetype = "solid") + 
  geom_label_repel(data = df.species,
                   mapping = aes(label= functional.group, x = MDS2, MDS3),
                   colour = "black",
                   alpha = 0.8,
                   nudge_x = 0,
                   nudge_y = 0,
                   size = label.size, 
                   box.padding = 0.5,
                   direction = "both",
                   max.overlaps = p.nr.functional.groups) +
  theme(legend.position="none") +
  annotate("text",
           -Inf, Inf,
           label = "Loadings",
           hjust = title.position.x,
           vjust = title.position.y,
           size = 5,
           fontface = "bold")


df.table = df.species[,c("ID", "functional.group")]

## Plot Legend
p4 = ggplot() + geom_point(data = subset(df.points, Scenario=="Current"),
                           aes(x= MDS2, y = MDS3, fill = Biogeo.realm,
                               col = Biogeo.realm),
                           alpha = 0, 
                           shape = 6,
                           size = 1.5) +
  geom_point(data = subset(df.points, Scenario=="Present-Natural"),
             aes(x= MDS2, y = MDS3, fill = Biogeo.realm,
                 col = Biogeo.realm),
             alpha = 0, 
             shape = 16,
             size = 2.5) +
  # annotate(geom = "table", x = min(df.points$MDS2), y = max(df.points$MDS3), label = list(df.table), 
  #          vjust = 1, hjust = 0, size = 3) +
  theme_void() +
  theme(legend.position="none")


# get_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}

legend = ggarrange(p4,p1_legend)

(plot.final.3 = ggarrange(p1 + theme(legend.position="none"),
                          p2 + theme(legend.position="none"),
                          p3 + theme(legend.position="none"),
                          legend,
                          ncol = 2,
                          nrow = 2,
                          labels = c("A", "B", "C")))


png(width=1000,height=1000, file = "3_Output/IAa7_Run_NMDS/IAa7_NMDS_Output_dim2_3.png")
plot.final.3
dev.off()


pdf(width=10,height=10, file = "3_Output/IAa7_Run_NMDS/IAa7_NMDS_Output_dim2_3.pdf")
plot.final.3
dev.off()





## – End Time ####
#–––––––––––––––––––––––––––––––––––––––––––––––––
end_time <- Sys.time()
end_time - start_time




    #~~~~~ SOME INFORMATION ON INTERPRETING NMDS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    #  SOURCE: http://stratigrafia.org/8370/lecturenotes/multidimensionalScaling.html
    #
    #   NMDS is not an eigenanalysis technique like principal components analysis or correspondence analysis. As a result, the axes cannot be interpreted such that axis 1 explains the greatest amount of variance, axis 2 explains the next greatest amount of variance, and so on. As a result, an NMS ordination can be rotated, inverted, or centered to any desired configuration.
    #   Benefits of NMDS: Rank-order (non-metric) approach well-suited for certain types of data (particularly counts of abundance). Uses specialized (and user-selected) distance metrics designed for data, such as Bray-Curtis (ecology) or Jaccard (biogeography), or any other.
    #   Next, a desired number of k dimensions is chosen for the ordination. The resulting ordination can be greatly sensitive to this number of chosen dimensions. For example, a k-dimensional ordination is not equivalent to the first k dimensions of a k+1-dimensional ordination.
    #   The ordination is sensitive to the number of dimensions that is chosen, so this choice must be made with care. Choosing too few dimensions will force multiple axes of variation to be expressed on a single ordination dimension. Choosing too many dimensions is no better in that it can cause a single source of variation to be expressed on more than one dimension. One way to choose an appropriate number of dimensions is perform ordinations of progressively higher numbers of dimensions. A scree diagram (stress versus number of dimensions) can then be plotted, on which one can identify the point beyond which additional dimensions do not substantially lower the stress value. A second criterion for the appropriate number of dimensions is the interpretability of the ordination, that is, whether the results make sense.
    #   The goal of NMDS is to represent the position of objects (e.g. communities) in multidimensional space as accurately as possible using a reduced number of dimensions that can be easily visualized (similar to PCA). The main aim is to recognize and interpret patterns and find gradients that represent the underlying geographical, ecological etc. gradients.
    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
