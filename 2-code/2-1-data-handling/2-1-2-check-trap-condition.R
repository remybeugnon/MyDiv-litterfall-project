#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-23 Check condition of traps - intact, broken, missing stones
# by Elisabeth Bönisch (elisabeth.boenisch@idiv.de)
# updated on 

#============================  Packages   ===============================

rm(list=ls())

library(readr)
library(readxl)
library(ggplot2)
library(tidyverse)

library(nlme)
library(lme4)
library(lmerTest)

library(devtools)
library(httr)
library(jsonlite)
library(XML)
install_github("BEXIS2/rBExIS", subdir = "rBExIS", dependencies=TRUE)
library(rBExIS)
load_all("rBExIS")
check("rBExIS")
bexis.options("base_url" = "https://mydivdata.idiv.de")
bexis.get.datasets()

setwd("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project")

#============================  Read the data   =======================================

d1_wide <- read.csv("C:/Users/eb64wupi/Documents/MyDiv/MyDiv-litterfall-project/1-data/2-1-data-handling/2-1-1-Full-data-wideformat-MyDiv-litter-dryweight.csv")



# throwing out outlier of 253.200??? or is it normal???
d1_wide$"Aesculus.hippocastanum"[d1_wide$"Aesculus.hippocastanum" == 253.200] <- NA

#### total ####
ggplot(d1_wide, aes(trap_condition, fill=trap_condition, color=trap_condition))+
  scale_fill_manual(values= c("#335C67","#E09F3E","#AD9BAA"))+ 
  scale_color_manual(values= c("black","black","black"))+
  geom_bar()+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))#

# per month
# change sequence of month
d1_wide$month1 = factor(d1_wide$month, levels=c("March","April","May","June","July","August","September","October","November","December","January","February"))

ggplot(d1_wide, aes(trap_condition, fill=trap_condition))+
  scale_fill_manual(values= c("#335C67","#E09F3E","#AD9BAA"))+ 
  scale_color_manual(values= c("black","black","black"))+
  geom_bar()+
  facet_grid(.~month1)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12,  angle = 90),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))

# per species
d1_long <- d1_wide %>%
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
                cont = "contaminants")%>% 
  pivot_longer(cols=c(8:18),
               names_to="species",
               values_to="dryweight") 

ggplot(d1_long, aes(trap_condition, fill=trap_condition))+
  geom_bar()+
  scale_fill_manual(values= c("#335C67","#E09F3E","#AD9BAA"))+ 
  scale_color_manual(values= c("black","black","black"))+
  facet_grid(.~species)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12,  angle = 90),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))

### end ###