#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-23
# by Elisabeth Bönisch (elisabeth.boenisch@idiv.de)
# updated on 

#============================  Clean the environment   ===============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)

library(nlme)
library(lme4)
library(lmerTest)

library(devtools)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
load_all("rBExIS")
bexis.options("base_url" = "https://mydivdata.idiv.de")
bexis.get.datasets()
# setwd("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project")
# .x = "2023-03-litter dry weights.xlsx"
df.all.wide = 
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
                  "Other (g)"), as.numeric) |> 
      dplyr::rename("trapID" = "Trap ID",
                    "plotID" = "Plot ID",
                    "other (g)" = "Other (g)") |>
      dplyr::select(-c(`Sampling date (YYYY-MM-DD)`)) |> 
      mutate(month = str_sub(.x, 1, 7))
  })

#============================  Read the data   =======================================



March23 <- read_excel("1-data/1-1-dry-weights/2023-03-litter dry weights.xlsx")
April23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-04-litter dry weights.xlsx")
May23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-05-litter dry weights.xlsx")
June23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-06-litter dry weights.xlsx")
July23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-07-litter dry weights.xlsx")
August23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-08-litter dry weights.xlsx")
September23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-09-litter dry weights.xlsx")
October23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-10-litter dry weights.xlsx")
November23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-11-litter dry weights.xlsx")
December23 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2023-12-litter dry weights.xlsx")
January24 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2024-01-litter dry weights.xlsx")
February24 <- read_excel("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/1-1-dry-weights/2024-02-litter dry weights.xlsx")

plotinfo <- rBExIS::bexis.GetDatasetById(id =43)

#============================  prepare data for joining   =======================================

# add months as variable
March23$month <- c("March")

March23 <- March23 %>%
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
              "Other (g)"), as.numeric)%>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID",
                "other (g)" = "Other (g)")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))

April23$month <- c("April")

April23<- April23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

April23 <- April23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))

May23$month <- c("May")

May23 <- May23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

May23 <- May23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))


June23$month <- c("June")

June23 <- June23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

June23 <- June23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))


July23$month <- c("July")

July23 <- July23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

July23 <- July23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`)) %>% 
  dplyr::filter(row_number() <= n()-1)


August23$month <- c("August")

August23 <- August23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

August23 <- August23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))


September23$month <- c("September")

September23 <- September23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

September23 <- September23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))

October23$month <- c("October")

October23 <- October23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

October23 <- October23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))

November23$month <- c("November")

November23 <- November23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

November23 <- November23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))

December23$month <- c("December")

December23 <- December23 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

December23 <- December23 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))

January24$month <- c("January")

January24 <- January24 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

January24 <- January24 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))

February24$month <- c("February")

February24 <- February24 %>%
  dplyr::mutate_at(c("Acer pseudoplatanus (g)",	
                     "Aesculus hippocastanum (g)",
                     "Betula pendula (g)",	
                     "Carpinus betulus (g)",
                     "Fagus sylvatica (g)",	
                     "Fraxinus excelsior (g)",	
                     "Prunus avium (g)",	
                     "Quercus petraea (g)", 
                     "Sorbus aucuparia (g)",	
                     "Tilia platyphyllos (g)", 
                     "other (g)"), as.numeric)

February24 <- February24 %>%
  dplyr::rename("trapID" = "Trap ID",
                "plotID" = "Plot ID")%>%
  dplyr::select(-c(`Sampling date (YYYY-MM-DD)`))


#============================  join the months  =======================================
# March + April
MarchApril <- full_join(March23, April23, by=c("plotID",
                                               "trapID",
                                               "Acer pseudoplatanus (g)",	
                                               "Aesculus hippocastanum (g)",
                                               "Betula pendula (g)",	
                                               "Carpinus betulus (g)",
                                               "Fagus sylvatica (g)",	
                                               "Fraxinus excelsior (g)",	
                                               "Prunus avium (g)",	
                                               "Quercus petraea (g)", 
                                               "Sorbus aucuparia (g)",	
                                               "Tilia platyphyllos (g)", 
                                               "other (g)",
                                               "trap condition",
                                               "month"))
# March + April + May
MAAPMA <- full_join(MarchApril, May23, by=c("plotID",
                                            "trapID",
                                            "Acer pseudoplatanus (g)",	
                                            "Aesculus hippocastanum (g)",
                                            "Betula pendula (g)",	
                                            "Carpinus betulus (g)",
                                            "Fagus sylvatica (g)",	
                                            "Fraxinus excelsior (g)",	
                                            "Prunus avium (g)",	
                                            "Quercus petraea (g)", 
                                            "Sorbus aucuparia (g)",	
                                            "Tilia platyphyllos (g)", 
                                            "other (g)",
                                            "trap condition",
                                            "month"))

# March + April + May + June
MAAPMAJ <- full_join(MAAPMA, June23, by=c("plotID",
                                          "trapID",
                                          "Acer pseudoplatanus (g)",	
                                          "Aesculus hippocastanum (g)",
                                          "Betula pendula (g)",	
                                          "Carpinus betulus (g)",
                                          "Fagus sylvatica (g)",	
                                          "Fraxinus excelsior (g)",	
                                          "Prunus avium (g)",	
                                          "Quercus petraea (g)", 
                                          "Sorbus aucuparia (g)",	
                                          "Tilia platyphyllos (g)", 
                                          "other (g)",
                                          "trap condition",
                                          "month"))

# March + April + May + June + July
MAAPMAJJ <- full_join(MAAPMAJ, July23, by=c("plotID",
                                            "trapID",
                                            "Acer pseudoplatanus (g)",	
                                            "Aesculus hippocastanum (g)",
                                            "Betula pendula (g)",	
                                            "Carpinus betulus (g)",
                                            "Fagus sylvatica (g)",	
                                            "Fraxinus excelsior (g)",	
                                            "Prunus avium (g)",	
                                            "Quercus petraea (g)", 
                                            "Sorbus aucuparia (g)",	
                                            "Tilia platyphyllos (g)", 
                                            "other (g)",
                                            "trap condition",
                                            "month"))

# March + April + May + June + July + August
MAAPMAJJA <- full_join(MAAPMAJJ, August23, by=c("plotID",
                                                "trapID",
                                                "Acer pseudoplatanus (g)",	
                                                "Aesculus hippocastanum (g)",
                                                "Betula pendula (g)",	
                                                "Carpinus betulus (g)",
                                                "Fagus sylvatica (g)",	
                                                "Fraxinus excelsior (g)",	
                                                "Prunus avium (g)",	
                                                "Quercus petraea (g)", 
                                                "Sorbus aucuparia (g)",	
                                                "Tilia platyphyllos (g)", 
                                                "other (g)",
                                                "trap condition",
                                                "month"))

# March + April + May + June + July + August + September
MAAPMAJJAS <- full_join(MAAPMAJJA, September23, by=c("plotID",
                                                     "trapID",
                                                     "Acer pseudoplatanus (g)",	
                                                     "Aesculus hippocastanum (g)",
                                                     "Betula pendula (g)",	
                                                     "Carpinus betulus (g)",
                                                     "Fagus sylvatica (g)",	
                                                     "Fraxinus excelsior (g)",	
                                                     "Prunus avium (g)",	
                                                     "Quercus petraea (g)", 
                                                     "Sorbus aucuparia (g)",	
                                                     "Tilia platyphyllos (g)", 
                                                     "other (g)",
                                                     "trap condition",
                                                     "month"))

# March + April + May + June + July + August + September + October
MAAPMAJJASO <- full_join(MAAPMAJJAS, October23, by=c("plotID",
                                                     "trapID",
                                                     "Acer pseudoplatanus (g)",	
                                                     "Aesculus hippocastanum (g)",
                                                     "Betula pendula (g)",	
                                                     "Carpinus betulus (g)",
                                                     "Fagus sylvatica (g)",	
                                                     "Fraxinus excelsior (g)",	
                                                     "Prunus avium (g)",	
                                                     "Quercus petraea (g)", 
                                                     "Sorbus aucuparia (g)",	
                                                     "Tilia platyphyllos (g)", 
                                                     "other (g)",
                                                     "trap condition",
                                                     "month"))

# March + April + May + June + July + August + September + October + November
MAAPMAJJASON <- full_join(MAAPMAJJASO, November23, by=c("plotID",
                                                        "trapID",
                                                        "Acer pseudoplatanus (g)",	
                                                        "Aesculus hippocastanum (g)",
                                                        "Betula pendula (g)",	
                                                        "Carpinus betulus (g)",
                                                        "Fagus sylvatica (g)",	
                                                        "Fraxinus excelsior (g)",	
                                                        "Prunus avium (g)",	
                                                        "Quercus petraea (g)", 
                                                        "Sorbus aucuparia (g)",	
                                                        "Tilia platyphyllos (g)", 
                                                        "other (g)",
                                                        "trap condition",
                                                        "month"))

# March + April + May + June + July + August + September + October + November + December
MAAPMAJJASOND <- full_join(MAAPMAJJASON, December23, by=c("plotID",
                                                          "trapID",
                                                          "Acer pseudoplatanus (g)",	
                                                          "Aesculus hippocastanum (g)",
                                                          "Betula pendula (g)",	
                                                          "Carpinus betulus (g)",
                                                          "Fagus sylvatica (g)",	
                                                          "Fraxinus excelsior (g)",	
                                                          "Prunus avium (g)",	
                                                          "Quercus petraea (g)", 
                                                          "Sorbus aucuparia (g)",	
                                                          "Tilia platyphyllos (g)", 
                                                          "other (g)",
                                                          "trap condition",
                                                          "month"))

# March + April + May + June + July + August + September + October + November + December + January
MAAPMAJJASONDJ <- full_join(MAAPMAJJASOND, January24, by=c("plotID",
                                                           "trapID",
                                                           "Acer pseudoplatanus (g)",	
                                                           "Aesculus hippocastanum (g)",
                                                           "Betula pendula (g)",	
                                                           "Carpinus betulus (g)",
                                                           "Fagus sylvatica (g)",	
                                                           "Fraxinus excelsior (g)",	
                                                           "Prunus avium (g)",	
                                                           "Quercus petraea (g)", 
                                                           "Sorbus aucuparia (g)",	
                                                           "Tilia platyphyllos (g)", 
                                                           "other (g)",
                                                           "trap condition",
                                                           "month"))

# March + April + May + June + July + August + September + October + November + December + January + February
MAAPMAJJASONDJF <- full_join(MAAPMAJJASONDJ, February24, by=c("plotID",
                                                              "trapID",
                                                              "Acer pseudoplatanus (g)",	
                                                              "Aesculus hippocastanum (g)",
                                                              "Betula pendula (g)",	
                                                              "Carpinus betulus (g)",
                                                              "Fagus sylvatica (g)",	
                                                              "Fraxinus excelsior (g)",	
                                                              "Prunus avium (g)",	
                                                              "Quercus petraea (g)", 
                                                              "Sorbus aucuparia (g)",	
                                                              "Tilia platyphyllos (g)", 
                                                              "other (g)",
                                                              "trap condition",
                                                              "month"))

# join with plot info (from MyDiv database)
litter_dryweights_fullyear_wide <- full_join(plotinfo, MAAPMAJJASONDJF, by=c("plotName"="plotID"))

# rename the species
litter_dryweights_fullyear_wide <- litter_dryweights_fullyear_wide %>%
  dplyr::rename("Acer pseudoplatanus" = "Acer pseudoplatanus (g)",	
                "Aesculus hippocastanum" = "Aesculus hippocastanum (g)",
                "Betula pendula" = "Betula pendula (g)",	
                "Carpinus betulus" = "Carpinus betulus (g)",
                "Fagus sylvatica" = "Fagus sylvatica (g)",	
                "Fraxinus excelsior" = "Fraxinus excelsior (g)",	
                "Prunus avium" = "Prunus avium (g)",	
                "Quercus petraea" = "Quercus petraea (g)", 
                "Sorbus aucuparia" = "Sorbus aucuparia (g)",	
                "Tilia platyphyllos" = "Tilia platyphyllos (g)", 
                "contaminants" = "other (g)")

#============================  add and change variable names/sequence   ======================================= 

litter_dryweights_fullyear_wide <- litter_dryweights_fullyear_wide %>%
  dplyr::rename("trap_condition"="trap condition")
litter_dryweights_fullyear_wide$div<-as.factor(litter_dryweights_fullyear_wide$tree_species_richness)
litter_dryweights_fullyear_wide$blk<-as.factor(litter_dryweights_fullyear_wide$block)
litter_dryweights_fullyear_wide$myc<-as.factor(litter_dryweights_fullyear_wide$mycorrhizal_type)
litter_dryweights_fullyear_wide$myc <- recode_factor(litter_dryweights_fullyear_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")
litter_dryweights_fullyear_wide$sr<-litter_dryweights_fullyear_wide$tree_species_richness
litter_dryweights_fullyear_wide$sr_myc<-paste(litter_dryweights_fullyear_wide$sr,litter_dryweights_fullyear_wide$myc,sep="_")#interaction terms
litter_dryweights_fullyear_wide$month1 = factor(litter_dryweights_fullyear_wide$month, levels=c("March","April","May","June","July","August","September","October","November","December","January","February"))


#============================  save full datatable  ======================================= 

write.csv(litter_dryweights_fullyear_wide, "C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/2-1-data-handling/2-1-1-Full-data-wideformat-MyDiv-litter-dryweight.csv", row.names = FALSE)

### end

