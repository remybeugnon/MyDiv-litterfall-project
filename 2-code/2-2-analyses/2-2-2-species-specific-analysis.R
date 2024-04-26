#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-23 Species specific litter fall data
# by Elisabeth Bönisch (elisabeth.boenisch@idiv.de)
# updated on ###

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

# convert to long format for plotting species 

d1_long <- d1_wide %>% 
  pivot_longer(cols=c(8:18),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich plots 
d1_long <- d1_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# change sequence of months
d1_long$month1 = factor(d1_long$month1, levels=c("March","April","May","June","July","August","September","October","November","December","January","February"))


# plot the species specific litter dryweight
ggplot(d1_long, aes(x=myc, y=dryweight_corr, color =species))+
  geom_boxplot()+
  geom_jitter(shape=21, size=2, alpha=0.4)+
  facet_grid(month1~div)+
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
        legend.box.background = element_rect(color="transparent"))


# regression lines

ggplot(d1_long, aes(x=sr, y=dryweight_corr, color =species))+
  geom_point(shape=21, size=2, alpha=0.4)+
  geom_smooth(method="lm", fullrange = FALSE, se = T, aes(colour= sr))+
  facet_grid(month~myc)+
  scale_x_continuous(breaks=c(1,2,4),
                     labels = c("1","2","4"))+
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
        legend.box.background = element_rect(color="transparent"))
