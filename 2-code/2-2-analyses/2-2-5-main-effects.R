#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-30
# main effects (tree diversity and mycorrhizal type)
# by Elisabeth Bönisch (elisabeth.boenisch@idiv.de)
# updated on...

#============================ Packages ===============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(lmerTest)

library(devtools)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
load_all("rBExIS")
bexis.options("base_url" = "https://mydivdata.idiv.de")
bexis.get.datasets()

#============================ Dataset ===============================

df.all.wide.info <- read_csv("1-data/2-1-data-handling/2-1-1-Full-data-wideformat-MyDiv-litter-dryweight.csv")

# 1) average across traps ####
df.all.wide.mean <- df.all.wide.info %>%
  dplyr::rename(Ac = "Acer pseudoplatanus",	
                Ae = "Aesculus hippocastanum",
                Be = "Betula pendula",	
                Ca = "Carpinus betulus",
                Fa = "Fagus sylvatica",	
                Fr = "Fraxinus excelsior",	
                Pr = "Prunus avium",	
                Qu = "Quercus petraea", 
                So = "Sorbus aucuparia",	
                Ti = "Tilia platyphyllos", 
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

# 2) long format ####
df.all.long <- df.all.wide.mean %>% 
  pivot_longer(cols=c(13:22),
               names_to="species",
               values_to="dryweight")


# 3) plot main interaction effects ####
df.annual.litter = 
  df.all.long |> 
  group_by(plotID, myc, tree_species_richness, block, composition) |>
  summarise(litter.prod = sum(dryweight, na.rm = T), 
            sd.litter.prod = sd(dryweight, na.rm = T))

ggplot(df.annual.litter, 
       aes(x=tree_species_richness, y=litter.prod, 
           color = myc, fill = myc))+
  geom_point(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
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

# 4) Model ####
mod.total.litterfall =
  lmerTest::lmer(litter.prod ~ myc * tree_species_richness + 
                   (1|block),
                 data = df.annual.litter)

# 5) Check the model quality ####
library(performance)
png("3-plots/2-2-5-Check-model-main-effects-2024-05-06.png")
performance::check_model(mod.total.litterfall)
dev.off()

# 6) Summary ####
summary(mod.total.litterfall)

# /) Anova (Type III SS) ####
anova(mod.total.litterfall)


### end ###

