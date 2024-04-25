#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-25
# diversity effects 
# by Elisabeth Bönisch (elisabeth.boenisch@idiv.de)
# updated on 

#============================ Packages ===============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)

library(devtools)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
load_all("rBExIS")
bexis.options("base_url" = "https://mydivdata.idiv.de")
bexis.get.datasets()


df.all.wide.info <- read.csv("2-1-1-Full-data-wideformat-MyDiv-litter-dryweight.csv")

df.all.long.info <- df.all.wide.info %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")
  
# correct for species richness
df.all.long.info <- df.all.long.info %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))       

# save as csv file 
write.csv(df.all.long.info, "1-data/2-1-data-handling/2-1-2-Full-data-longformat_corrected-MyDiv-litter-dryweight.csv", row.names = FALSE)

