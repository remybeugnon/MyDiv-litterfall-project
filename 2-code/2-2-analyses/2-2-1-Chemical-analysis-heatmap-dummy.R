#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-05-23
# chemical composition
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

#============================ Dataset ===============================

df.all.wide.info <- read.csv("1-data/2-1-data-handling/2-1-1-Full-data-wideformat-MyDiv-litter-dryweight-m2.csv")

# 1) average across traps ####
df.all.wide.mean <- df.all.wide.info %>%
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


# 2) long format ####
df.all.long <- df.all.wide.mean %>% 
  pivot_longer(cols=c(13:22),
               names_to="species",
               values_to="dryweight")


# 3) monthly effect - sum per plot and month ####
df.monthly.litter = df.all.long |> 
  group_by(block, plotID, sr, div, myc, month1) |> 
  summarise(litter.prod = sum(dryweight, na.rm = T))

df.monthly.litter$month2 <- dplyr::recode_factor(df.monthly.litter$month1, 
                                                 "March"="Mar", 
                                                 "April"="Apr",
                                                 "May"="May",
                                                 "June"="Jun",
                                                 "July"="Jul",
                                                 "August"="Aug",
                                                 "September"="Sep",
                                                 "October"="Oct",
                                                 "November"="Nov",
                                                 "December"="Dec",
                                                 "January"="Jan",
                                                 "February"="Feb")

df.monthly.litter$month1 = factor(df.monthly.litter$month2, 
                                  levels = c(month.abb[3:12], month.abb[1:2]))
# 
# df = 
#   tibble(
#     C = runif(120),
#     N = runif(120),
#     month = factor(month.name) |> rep(10)
#   ) |>
#   pivot_longer(cols = 1:2)


ggplot(data = df, 
       aes(x = month, y = name, fill = value)) + 
  geom_raster()
