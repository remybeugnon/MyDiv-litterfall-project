#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-05-14
# litter dryweight NMDS on species level over time
# by Elisabeth Bönisch (elisabeth.boenisch@idiv.de)
# updated on 

#============================ Packages ===============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(vegan)

#============================ Dataset ===============================

df.all.wide.info <- read.csv("2-1-1-Full-data-wideformat-MyDiv-litter-dryweight.csv")

# 1) average across traps ####
df.all.wide.mean <- df.all.wide.info %>%
  dplyr::rename(Ac = "Acer.pseudoplatanus",	
                Ae = "Aesculus.hippocastanum",
                Be = "Betula.pendula",	
                Ca = "Carpinus.betulus",
                Fa = "Fagus.sylvatica",	
                Fr = "Fraxinus.excelsior",	
                Pr = "Prunus.avium",	
                Qu = "Quercus.petraea", 
                So = "Sorbus.aucuparia",	
                Ti = "Tilia.platyphyllos", 
                cont = "contaminants") %>%
  dplyr::group_by(plotID, plotName, tree_species_richness, mycorrhizal_type, myc, sr, div, block, blk, month1, month, composition) %>% # removed trap and month
  dplyr::summarise(Ac_mean = mean(Ac, na.rm = TRUE),	
                   Ae_mean = mean(Ae, na.rm = TRUE),
                   Be_mean = mean(Be, na.rm = TRUE),	
                   Ca_mean = mean(Ca, na.rm = TRUE),
                   Fa_mean = mean(Fa, na.rm = TRUE),	
                   Fr_mean = mean(Fr, na.rm = TRUE),	
                   Pr_mean = mean(Pr, na.rm = TRUE),	
                   Qu_mean = mean(Qu, na.rm = TRUE), 
                   So_mean = mean(So, na.rm = TRUE),	
                   Ti_mean = mean(Ti, na.rm = TRUE), 
                   cont_mean = mean(cont, na.rm = TRUE))

df.all.wide.mean$div<-as.factor(df.all.wide.mean$tree_species_richness)
df.all.wide.mean$blk<-as.factor(df.all.wide.mean$block)
df.all.wide.mean$myc<-as.factor(df.all.wide.mean$mycorrhizal_type)
df.all.wide.mean$sr<-df.all.wide.mean$tree_species_richness
df.all.wide.mean$sr_myc<-paste(df.all.wide.mean$sr,df.all.wide.mean$myc,sep="_")#interaction terms
df.all.wide.mean$myc <- recode_factor(df.all.wide.mean$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")




# 2) plot by species matrix ####
# To run an ordination you will need a data-frame consisting of plot by species (or trait) matrix 

plot.by.species <- df.all.wide.mean |>
  dplyr::group_by(plotID, plotName, tree_species_richness, mycorrhizal_type, myc, sr, div, block, blk, month1, month, composition, Ac_mean, Ae_mean, Be_mean, Ca_mean, Fa_mean, Fr_mean, Pr_mean, Qu_mean, So_mean, Ti_mean) |>
  ungroup()|>
  dplyr::select(c("Ac_mean":"Ti_mean"))

# 3) groups dataframe ####
# AND a “groups” data-frame which should consist of plots with a coding variable for what group each plot belongs to, 
# this will be used for plotting the ordination.

plot.info <- df.all.wide.mean |>
  dplyr::group_by(plotID, plotName, tree_species_richness, mycorrhizal_type, myc, sr, div, block, blk, month1, month, composition, Ac_mean, Ae_mean, Be_mean, Ca_mean, Fa_mean, Fr_mean, Pr_mean, Qu_mean, So_mean, Ti_mean) |>
  ungroup()|>
  dplyr::select(c("plotID":"composition"))

# 4) Converting Absolute Abundance to Relative Abundance ####

tree.species.rel <-         
  decostand(plot.by.species, method = "total", na.rm = T)

# 5) Calculating your distance matrix ####
# Removing rows having all zeros
tree.species.rel2 <- tree.species.rel[rowSums(tree.species.rel[])>0,]

tree.species.distmat <- 
  vegdist(tree.species.rel2, method = "bray", na.rm = T)    # zero values are problem

# 6) Creating easy to view matrix and writing .csv ####
tree.species.distmat2 <- 
  as.matrix(tree.species.distmat, labels = T)
write.csv(tree.species.distmat, "1-data/2-2-9-NDMS-tree-species-distance-matrix.csv")


# 7) Running NMDS in vegan (metaMDS) ####                  # zero values & NAs are problem
tree.species_NMS <-
  metaMDS(tree.species.distmat2,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)
