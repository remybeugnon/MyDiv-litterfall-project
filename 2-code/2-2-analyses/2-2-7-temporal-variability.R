#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-05-11
# temporal variability of litter fall
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


# 3) temporal variability - sum per plot ####

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

#ggsave("3-plots/2-2-7-Figure-litterfall_total-per-trap-2024-05-11.jpeg", 
       total.tr, 
       height=20,
       width=28, 
       unit="cm", 
       dpi=2000) 

total.trap <- ggplot(df.1, aes(x=sr, y=litterfall_sum, color = myc, fill = myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5) +
  geom_smooth(method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Leaf litter dryweight total"~~(g/m2)), 
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

#ggsave("3-plots/2-2-7-Figure-litterfall_total-alltraps-2024-05-11.jpeg", 
       total.trap, 
       height=20,
       width=28, 
       unit="cm", 
       dpi=2000) 

### remove plot 50 ####
df.1.no50 = df.1[!(df.1$plotID == 50),]

total.trap.no50 <- ggplot(df.1.no50, aes(x=sr, y=litterfall_sum, color = myc, fill = myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5) +
  geom_smooth(method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Leaf litter dryweight total"~~(g/m2)), 
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

ggsave("3-plots/2-2-7-Figure-litterfall_total-alltraps-noplot50-2024-05-11.jpeg", 
total.trap.no50, 
height=20,
width=28, 
unit="cm", 
dpi=2000) 



# 5) temporal variability - sd of total litterfall   ####
#  Yearly spatial stability
df.2 = df.1 |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_sd = mean(litterfall_sum, na.rm=T)/sd(litterfall_sum, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID) |> 
  summarise(litterfall_sd = mean(litterfall_sd, na.rm=T))

#  Monthly spatial stability
df.2 = df.1 |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_sd = mean(litterfall_sum, na.rm=T)/sd(litterfall_sum, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2])))

# Yearly temporal stability 
df.2 = df.1 |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_mean = mean(litterfall_sum, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID) |> 
  summarise(litterfall_sd = mean(litterfall_mean, na.rm=T)/sd(litterfall_mean, na.rm=T))

# 6) plot of temporal variability - sd of total litterfall ####

temp.stab <-
ggplot()+
  geom_point(data = df.2, 
             aes(x=sr, y=litterfall_sd, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.2, 
              aes(x=sr, y=litterfall_sd, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Temporal stability of leaf litter production"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr'), 
            aes(x=1, y=305,
                label = sign,
                fontface = "bold",
                hjust = 0), 
            color = 'black') + 
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'myc'), 
            aes(x=1, y=295,
                label = sign,
                fontface = "bold",
                hjust = 0), 
            color = 'black') + 
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr:myc'), 
            aes(x=1, y=285,
                label = sign,
                fontface = "bold",
                hjust = 0), 
            color = 'black') + 
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

ggsave("3-plots/2-2-7-Figure-litterfalll-temp-stability-sig-2024-05-13.jpeg", 
       temp.stab, 
       height=20,
       width=28, 
       unit="cm", 
       dpi=2000) 

### remove plot 50 ####
df.2.no50 = df.2[!(df.2$plotID == 50),]

temp.stab.no50 <-ggplot(df.2.no50, aes(x=sr, y=litterfall_sd, color = myc, fill = myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5) +
  geom_smooth(method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Temporal stability of leaf litter production"), 
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

ggsave("3-plots/2-2-7-Figure-litterfalll-temp-stability-noplot50-2024-05-13.jpeg", 
temp.stab.no50, 
height=20,
width=28, 
unit="cm", 
dpi=2000) 


#7) Model ####

# with correlation structure
#library(nlme)
mod.temp.litter =
  lme(litterfall_sd ~ month1 * sr * myc,
      random = ~1|block,
      data = df.2,
      correlation=corCAR1())

# 6) Check the model quality ####
#library(performance)
png("3-plots/2-2-7-Check-model-temporal-stability-2024-05-13.png", 
    width=1000, height=1000)
performance::check_model(mod.temp.litter)
dev.off()

# 7) Summary ####
summary(mod.temp.litter)

# 8) Anova (Type III SS) ####
anova(mod.temp.litter)

# 9) Check months individually ####
M = map_df( .x = unique(df.2$month1),
            .f = ~ {
              mod = 
                lme(litterfall_sd ~ sr * myc, 
                    random= ~1|block,
                    data= df.2 |>
                      filter(month1 == .x)) |> 
                anova() |> 
                data.frame() |> 
                mutate(month = .x)
              mod$exp = rownames(mod)
              rownames(mod) <- NULL
              mod |> 
                select(month, explanatory = exp, 
                       numDF, denDF, F.value, p.value)
            }) |> 
  mutate(sign = if_else(p.value < 0.001, '***', 
                        if_else(p.value < 0.01, "**", 
                                if_else(p.value < 0.05, '*',
                                        if_else(p.value<0.1, '.', 'n.s.')))))
M

M$sign[M$explanatory == 'sr' & M$month == 'March'] = 
  paste0("sr = ",M$sign[M$explanatory == 'sr' & M$month == 'March'])
M$sign[M$explanatory == 'myc' & M$month == 'March'] = 
  paste0("myc = ",M$sign[M$explanatory == 'myc' & M$month == 'March'])
M$sign[M$explanatory == 'sr:myc' & M$month == 'March'] = 
  paste0("sr:myc = ",M$sign[M$explanatory == 'sr:myc' & M$month == 'March'])


