#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-23 Monthly analysis of litter dryweights
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


# extract the months


### ----- March ----- 
March_wide <- d1_wide %>%
  filter(month == "March")

March_wide$div<-as.factor(March_wide$tree_species_richness)
March_wide$blk<-as.factor(March_wide$block)
March_wide$myc<-as.factor(March_wide$mycorrhizal_type)
March_wide$sr<-March_wide$tree_species_richness
March_wide$sr_myc<-paste(March_wide$sr,March_wide$myc,sep="_")#interaction terms
March_wide$myc <- recode_factor(March_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
March_long <- March_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
March_long <- March_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

March <- ggplot(March_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(.~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

March

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-March.jpeg", 
March, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 

### ----- April ----- 
April_wide <- d1_wide %>%
  filter(month == "April")

April_wide$div<-as.factor(April_wide$tree_species_richness)
April_wide$blk<-as.factor(April_wide$block)
April_wide$myc<-as.factor(April_wide$mycorrhizal_type)
April_wide$sr<-April_wide$tree_species_richness
April_wide$sr_myc<-paste(April_wide$sr,April_wide$myc,sep="_")#interaction terms
April_wide$myc <- recode_factor(April_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
April_long <- April_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
April_long <- April_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

April <- ggplot(April_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(.~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

April

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-April.jpeg", 
April, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### ----- May ----- 
May_wide <- d1_wide %>%
  filter(month == "May")

May_wide$div<-as.factor(May_wide$tree_species_richness)
May_wide$blk<-as.factor(May_wide$block)
May_wide$myc<-as.factor(May_wide$mycorrhizal_type)
May_wide$sr<-May_wide$tree_species_richness
May_wide$sr_myc<-paste(May_wide$sr,May_wide$myc,sep="_")#interaction terms
May_wide$myc <- recode_factor(May_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
May_long <- May_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
May_long <- May_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

May <- ggplot(May_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(m.~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

May

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-May.jpeg", 
May, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### ----- June ----- 
June_wide <- d1_wide %>%
  filter(month == "June")

June_wide$div<-as.factor(June_wide$tree_species_richness)
June_wide$blk<-as.factor(June_wide$block)
June_wide$myc<-as.factor(June_wide$mycorrhizal_type)
June_wide$sr<-June_wide$tree_species_richness
June_wide$sr_myc<-paste(June_wide$sr,June_wide$myc,sep="_")#interaction terms
June_wide$myc <- recode_factor(June_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
June_long <- June_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
June_long <- June_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

June <- ggplot(June_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(.~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

June

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-June.jpeg", 
June, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 

### ----- July ----- 
July_wide <- d1_wide %>%
  filter(month == "July")

July_wide$div<-as.factor(July_wide$tree_species_richness)
July_wide$blk<-as.factor(July_wide$block)
July_wide$myc<-as.factor(July_wide$mycorrhizal_type)
July_wide$sr<-July_wide$tree_species_richness
July_wide$sr_myc<-paste(July_wide$sr,July_wide$myc,sep="_")#interaction terms
July_wide$myc <- recode_factor(July_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
July_long <- July_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
July_long <- July_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

July <- ggplot(July_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(.~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

July

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-July.jpeg", 
July, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### ----- August ----- 
Aug_wide <- d1_wide %>%
  filter(month == "August")

Aug_wide$div<-as.factor(Aug_wide$tree_species_richness)
Aug_wide$blk<-as.factor(Aug_wide$block)
Aug_wide$myc<-as.factor(Aug_wide$mycorrhizal_type)
Aug_wide$sr<-Aug_wide$tree_species_richness
Aug_wide$sr_myc<-paste(Aug_wide$sr,Aug_wide$myc,sep="_")#interaction terms
Aug_wide$myc <- recode_factor(Aug_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
Aug_long <- Aug_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
Aug_long <- Aug_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

Aug <- ggplot(Aug_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(.~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

Aug

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-Aug.jpeg", 
Aug, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 



### ----- September ----- 
Sep_wide <- d1_wide %>%
  filter(month == "September")

Sep_wide$div<-as.factor(Sep_wide$tree_species_richness)
Sep_wide$blk<-as.factor(Sep_wide$block)
Sep_wide$myc<-as.factor(Sep_wide$mycorrhizal_type)
Sep_wide$sr<-Sep_wide$tree_species_richness
Sep_wide$sr_myc<-paste(Sep_wide$sr,Sep_wide$myc,sep="_")#interaction terms
Sep_wide$myc <- recode_factor(Sep_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
Sep_long <- Sep_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
Sep_long <- Sep_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

Sep <- ggplot(Sep_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(month~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

Sep

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-Sep.jpeg", 
Sep, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### ----- October ----- 
Oct_wide <- d1_wide %>%
  filter(month == "October")

Oct_wide$div<-as.factor(Oct_wide$tree_species_richness)
Oct_wide$blk<-as.factor(Oct_wide$block)
Oct_wide$myc<-as.factor(Oct_wide$mycorrhizal_type)
Oct_wide$sr<-Oct_wide$tree_species_richness
Oct_wide$sr_myc<-paste(Oct_wide$sr,Oct_wide$myc,sep="_")#interaction terms
Oct_wide$myc <- recode_factor(Oct_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
Oct_long <- Oct_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
Oct_long <- Oct_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

Oct <- ggplot(Oct_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(month~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

Oct

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-Oct.jpeg", 
Oct, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### ----- November ----- 
Nov_wide <- d1_wide %>%
  filter(month == "November")

Nov_wide$div<-as.factor(Nov_wide$tree_species_richness)
Nov_wide$blk<-as.factor(Nov_wide$block)
Nov_wide$myc<-as.factor(Nov_wide$mycorrhizal_type)
Nov_wide$sr<-Nov_wide$tree_species_richness
Nov_wide$sr_myc<-paste(Nov_wide$sr,Nov_wide$myc,sep="_")#interaction terms
Nov_wide$myc <- recode_factor(Nov_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
Nov_long <- Nov_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
Nov_long <- Nov_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

Nov <- ggplot(Nov_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(month~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

Nov

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-Nov.jpeg", 
Nov, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### ----- December ----- 
Dec_wide <- d1_wide %>%
  filter(month == "December")

Dec_wide$div<-as.factor(Dec_wide$tree_species_richness)
Dec_wide$blk<-as.factor(Dec_wide$block)
Dec_wide$myc<-as.factor(Dec_wide$mycorrhizal_type)
Dec_wide$sr<-Dec_wide$tree_species_richness
Dec_wide$sr_myc<-paste(Dec_wide$sr,Dec_wide$myc,sep="_")#interaction terms
Dec_wide$myc <- recode_factor(Dec_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
Dec_long <- Dec_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
Dec_long <- Dec_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

Dec <- ggplot(Dec_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(month~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

Dec

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-Dec.jpeg", 
Dec, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### ----- January ----- 
Jan_wide <- d1_wide %>%
  filter(month == "January")

Jan_wide$div<-as.factor(Jan_wide$tree_species_richness)
Jan_wide$blk<-as.factor(Jan_wide$block)
Jan_wide$myc<-as.factor(Jan_wide$mycorrhizal_type)
Jan_wide$sr<-Jan_wide$tree_species_richness
Jan_wide$sr_myc<-paste(Jan_wide$sr,Jan_wide$myc,sep="_")#interaction terms
Jan_wide$myc <- recode_factor(Jan_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
Jan_long <- Jan_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
Jan_long <- Jan_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

Jan <- ggplot(Jan_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  facet_grid(month~div)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

Jan

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-Jan.jpeg", 
Jan, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### ----- February ----- 
Feb_wide <- d1_wide %>%
  filter(month == "February")

Feb_wide$div<-as.factor(Feb_wide$tree_species_richness)
Feb_wide$blk<-as.factor(Feb_wide$block)
Feb_wide$myc<-as.factor(Feb_wide$mycorrhizal_type)
Feb_wide$sr<-Feb_wide$tree_species_richness
Feb_wide$sr_myc<-paste(Feb_wide$sr,Feb_wide$myc,sep="_")#interaction terms
Feb_wide$myc <- recode_factor(Feb_wide$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# convert to long format for plotting
Feb_long <- Feb_wide %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")

# correct for lower number of individuals per species in species rich communities
Feb_long <- Feb_long %>%
  dplyr::mutate(dryweight_corr = case_when(tree_species_richness == 1 ~ dryweight*1,
                                           tree_species_richness == 2 ~ dryweight*2,
                                           tree_species_richness == 4 ~ dryweight*4))

# plot the species specific litter dryweight

Feb <- ggplot(Feb_long, aes(x=myc, y=dryweight_corr, color =species, fill =species))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.4)+
  labs(y=bquote("Leaf litter dryweight"~~(g)),   ### three traps, not averaged
       x = "Tree species richness")+
  # scale_x_continuous(breaks=c(1,2,4))+
  #scale_y_continuous(limits=c(15,35))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate(geom="text", x=3.5, y=32, label="italic(p) = 0.538",
  #          color="red", size = 7)+
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

Feb

#ggsave("Fig-2-1-4-MyDiv-litter-dryweight-specieslevel-boxplots-Feb.jpeg", 
Feb, 
height=12,
width=22, 
unit="cm", 
dpi=1000) 


### end ###

