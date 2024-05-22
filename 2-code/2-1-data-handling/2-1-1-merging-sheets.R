#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# merging sheets
# 2024-04-25
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

plotinfo <- rBExIS::bexis.GetDatasetById(id =43)

#setwd("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project")
.x = "2023-03-litter dry weights.xlsx"
df.all.wide.m2 = 
  list.files('1-data/1-1-dry-weights') |> 
  map_df(.f = ~{
    d = read_excel(paste0("1-data/1-1-dry-weights/",.x)) |>
      mutate_at(c("Acer pseudoplatanus (g)",	
                  "Aesculus hippocastanum (g)",
                  "Betula pendula (g)",	
                  "Carpinus betulus (g)",
                  "Fagus sylvatica (g)",	
                  "Fraxinus excelsior (g)",	
                  "Prunus avium (g)",	
                  "Quercus petraea (g)", 
                  "Sorbus aucuparia (g)",	
                  "Tilia platyphyllos (g)", 
                  "contaminants (g)"), as.numeric) |> 
      dplyr::rename("trapID" = "Trap ID",
                    "plotID" = "Plot ID",
                    "Ac" = "Acer pseudoplatanus (g)",	
                    "Ae" = "Aesculus hippocastanum (g)",
                    "Be" = "Betula pendula (g)",	
                    "Ca" = "Carpinus betulus (g)",
                    "Fa" = "Fagus sylvatica (g)",	
                    "Fr" = "Fraxinus excelsior (g)",	
                    "Pr" = "Prunus avium (g)",	
                    "Qu" = "Quercus petraea (g)", 
                    "So" = "Sorbus aucuparia (g)",	
                    "Ti" = "Tilia platyphyllos (g)", 
                    "cont" = "contaminants (g)",
                    "trap_condition"="trap condition") |>
      dplyr::mutate(Ac_m2 = (Ac/0.5625),	
                    Ae_m2 = (Ae/0.5625),
                    Be_m2 = (Be/0.5625),
                    Ca_m2 = (Ca/0.5625),
                    Fa_m2 = (Fa/0.5625), 
                    Fr_m2 = (Fr/0.5625),
                    Pr_m2 = (Pr/0.5625),
                    Qu_m2 = (Qu/0.5625), 
                    So_m2 = (So/0.5625),
                    Ti_m2 = (Ti/0.5625), 
                    cont_m2 = (cont/0.5625)) |>
      dplyr::select(-c(`Sampling date (YYYY-MM-DD)`)) |> 
      mutate(month = str_sub(.x, 1, 7))
  })

  
# remove row 1201

df.all.wide.m2 <- df.all.wide.m2[-c(1201), ]

df.all.wide.info <- full_join(plotinfo, df.all.wide.m2, by=c("plotName"="plotID"))


#============================  add and change variable names/sequence   ======================================= 

df.all.wide.info$div<-as.factor(df.all.wide.info$tree_species_richness)
df.all.wide.info$blk<-as.factor(df.all.wide.info$block)
df.all.wide.info$myc<-as.factor(df.all.wide.info$mycorrhizal_type)
df.all.wide.info$myc<-recode_factor(df.all.wide.info$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")
df.all.wide.info$sr<-df.all.wide.info$tree_species_richness
df.all.wide.info$sr_myc<-paste(df.all.wide.info$sr,df.all.wide.info$myc,sep="_")#interaction terms

df.all.wide.info$month1<-as.factor(df.all.wide.info$month)
df.all.wide.info$month1<-recode_factor(df.all.wide.info$month1, "2023-03"="March","2023-04"="April","2023-05"="May","2023-06"="June","2023-07"="July","2023-08"="August","2023-09"="September","2023-10"="October","2023-11"="November","2023-12"="December","2024-01"="January","2024-02"="February" )


#============================  save full datatable  ======================================= 

write.csv(df.all.wide.info, "1-data/2-1-data-handling/2-1-1-Full-data-wideformat-MyDiv-litter-dryweight-m2.csv", row.names = FALSE)

### end

