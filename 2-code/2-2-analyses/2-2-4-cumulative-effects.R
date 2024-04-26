#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-26
# cumulative sum of biomass 
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

df.all.long = read_csv('1-data/2-1-data-handling/2-1-2-Full-data-longformat_corrected-MyDiv-litter-dryweight.csv')

df = df.all.long |> 
  group_by(block, plotID, div, myc, month1) |> 
  summarise(litterfall = sum(dryweight, na.rm = T))

df$month1 = factor(df$month1, levels = 
                     c( month.name[3:12],  month.name[1:2]))
df.1 = 
  df |> 
  group_by(block, plotID, div, myc) |> 
  arrange(month1) |> 
  mutate(cs = cumsum(litterfall))


ggplot(data = df.1, aes(x = as.numeric(month1), y = cs, color = factor(myc), fill = factor(myc))) + 
  geom_smooth(alpha=0.2)+
  geom_jitter(shape = 21, alpha=0.5)+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  # scale_x_discrete(df.1$month1, 
  #                    breaks = c (3,4,5,6,7,8,9,10,11,12,1,2),
  #                    labels = c (3,4,5,6,7,8,9,10,11,12,1,2))+
  scale_fill_manual(values= c("#71b540","#febf00","#4c8ecb"),
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#febf00","#4c8ecb"), 
                     name = "Mycorrhizal type")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
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


# remove extreme plot 50 (Aesculus monoculture)
df.1.no50 = df.1[!(df.1$plotID == 50),]


ggplot(data = df.1.no50, aes(x = as.numeric(month1), y = cs, color = factor(div), fill = factor(div))) + 
  geom_smooth(alpha=0.2)+
  geom_jitter(shape = 21, alpha=0.5)+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  scale_x_continuous(breaks = c (3,4,5,6,7,8,9,10,11,12,1,2),
                     labels = c (3,4,5,6,7,8,9,10,11,12,1,2))+
  scale_fill_manual(values= c("#335C67","#E09F3E","#f68080"),
                    guide="none")+ 
  scale_color_manual(values= c("#335C67","#E09F3E","#f68080"), 
                     name = "Tree species richness")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
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


####  bar plot  ####

### myccorhizal type
ggplot(data = df.1, aes(x = as.numeric(month1), y = cs, color = factor(myc), fill = factor(myc))) + 
  geom_bar(stat="identity")+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  scale_x_continuous(breaks = c (3,4,5,6,7,8,9,10,11,12,1,2),
                     labels = c (3,4,5,6,7,8,9,10,11,12,1,2))+
  scale_fill_manual(values= c("#71b540","#febf00","#4c8ecb"),
                    name = "Mycorrhizal type")+ 
  scale_color_manual(values = c("#71b540","#febf00","#4c8ecb"),
                     guide ="none")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
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

### tree species richness
ggplot(data = df.1, aes(x = as.numeric(month1), y = cs, color = factor(div), fill = factor(div))) + 
  geom_bar(stat="identity")+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  scale_x_continuous(breaks = c (3,4,5,6,7,8,9,10,11,12,1,2),
                     labels = c (3,4,5,6,7,8,9,10,11,12,1,2))+
  scale_fill_manual(values= c("#335C67","#E09F3E","#f68080"),
                    name = "Tree species richness")+ 
  scale_color_manual(values= c("#335C67","#E09F3E","#f68080"), 
                     guide="none")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
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
