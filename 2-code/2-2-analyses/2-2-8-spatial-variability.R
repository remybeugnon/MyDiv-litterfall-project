#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-05-13
# spatial variability of litter fall
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

df.all.wide.info <- read.csv("2-1-1-Full-data-wideformat-MyDiv-litter-dryweight.csv")

# 1) NO! average across traps ####
df.all.wide.trap <- df.all.wide.info %>%
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
                cont = "contaminants")

df.all.wide.trap$div<-as.factor(df.all.wide.trap$tree_species_richness)
df.all.wide.trap$blk<-as.factor(df.all.wide.trap$block)
df.all.wide.trap$myc<-as.factor(df.all.wide.trap$mycorrhizal_type)
df.all.wide.trap$sr<-df.all.wide.trap$tree_species_richness
df.all.wide.trap$sr_myc<-paste(df.all.wide.trap$sr,df.all.wide.trap$myc,sep="_")#interaction terms
df.all.wide.trap$myc <- recode_factor(df.all.wide.trap$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")


# 2) long format ####
df.all.long <- df.all.wide.trap %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="dryweight")


# 3) spatial variability - sum per plot ####

df.1 = df.all.long |> 
  group_by(block, sr, div, myc, plotID, trapID, month1) |> 
  summarise(litterfall_sum = sum(dryweight, na.rm = TRUE)) |>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2])))

# 4) plot of total litter per trap ####

total.tr <- ggplot(df.1, 
                   aes(x=sr, y=litterfall_sum, 
                       color = myc, fill = myc))+
  geom_point(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  facet_grid(.~trapID)+
  labs(y=bquote("Leaf litter dryweight"~~(g/m2)), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
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

#ggsave("3-plots/2-2-8-Figure-litterfall_spatial-variability-2024-05-13.jpeg", 
total.tr, 
height=20,
width=28, 
unit="cm", 
dpi=2000) 
