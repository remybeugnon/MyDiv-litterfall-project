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

df.all.wide.info <- read.csv("2-1-1-Full-data-wideformat-MyDiv-litter-dryweight.csv")


# 1) average across traps ####
df.all.wide.mean <- df.all.wide.info %>%
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


# 3) sum across species 
df.cum.litter = df.all.long |> 
  group_by(block, plotID, div, sr, myc, month1) |> 
  summarise(litterfall = sum(dryweight, na.rm = T))

df.cum.litter$month1 = factor(df.cum.litter$month1, levels = 
                     c( month.name[3:12],  month.name[1:2]))
df.cum.litter.1 = 
  df.cum.litter |> 
  group_by(block, plotID, div, sr, myc) |> 
  arrange(month1) |> 
  mutate(cs = cumsum(litterfall))


# 4) plot cumulative effects ####
ggplot(data =df.cum.litter.1, aes(x = as.numeric(month1), y = cs, color = factor(myc), fill = factor(myc))) + 
  geom_smooth(alpha=0.2)+
  geom_jitter(shape = 21, alpha=0.5)+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  scale_x_discrete(breaks = 1:12,
                   labels = c (3,4,5,6,7,8,9,10,11,12,1,2))+
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




ggplot()+
  geom_point(data = df.cum.litter.1, 
             aes(x=sr, y=cs, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.cum.litter.1, 
              aes(x=sr, y=cs, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Cumulative sum - Leaf litter dryweight"~(g/m2)), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_y_continuous(limits = c(0, 130))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr'), 
            aes(x=2, y=130,
                label = sign,
                fontface = "bold",
                hjust = 0), 
            color = 'black') + 
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'myc'), 
            aes(x=2, y=125,
                label = sign,
                fontface = "bold",
                hjust = 0), 
            color = 'black') + 
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr:myc'), 
            aes(x=2, y=120,
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

# 4) Model ####

# with correlation structure
#library(nlme)
mod.cum.litter =
  lme(cs ~ month1 * sr * myc,
      random = ~1|block,
      data = df.cum.litter.1,
      correlation=corCAR1())

# 5) Check the model quality ####
#library(performance)
png("3-plots/2-2-4-Check-model-cumulative-sum-effects-2024-05-07.png")
performance::check_model(mod.cum.litter)
dev.off()

# 6) Summary ####
summary(mod.cum.litter)

# 7) Anova (Type III SS) ####
anova(mod.cum.litter)

# 8) Check months individually
M = map_df( .x = unique(df.cum.litter.1$month1),
            .f = ~ {
              mod = 
                lme(cs ~ sr * myc, 
                    random= ~1|block,
                    data= df.cum.litter.1 |>
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

### end ###






ggplot(data = df.1, aes(x = as.numeric(month1), y = cs, color = factor(myc), fill = factor(myc))) + 
  geom_smooth(alpha=0.2)+
  geom_jitter(shape = 21, alpha=0.5)+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  scale_x_discrete(breaks = 1:12,
                  labels = c (3,4,5,6,7,8,9,10,11,12,1,2))+
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


# test: remove extreme plot 50 (Aesculus monoculture)
df.1.no50 = df.1[!(df.1$plotID == 50),]


ggplot(data = df.1, aes(x = as.numeric(month1), y = cs, color = factor(div), fill = factor(div))) + 
  # geom_smooth(alpha=0.2)+
  geom_jitter(data = df.1, aes(x = as.numeric(month1), y = cs, color = factor(div), fill = factor(div)),
              shape = 21, alpha=0.5)+
  geom_line(data = df.1 |> 
              group_by(month1, div) |>
              summarise(cs = mean(cs)), 
            aes(x = as.numeric(month1), 
                y = cs, color = factor(div), fill = factor(div)),
              linewidth = 2)+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  scale_x_continuous(breaks = 1:12,
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
ggplot(data = df.1.no50, aes(x = as.numeric(month1), y = cs, color = factor(myc), fill = factor(myc))) + 
  geom_bar(position = "dodge", stat = "summary", fun = "mean"
  )+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  scale_x_continuous(breaks = 1:12,
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
ggplot(data = df.1.no50, aes(x = as.numeric(month1), y = cs, color = factor(div), fill = factor(div))) + 
  geom_bar(position = "dodge", stat = "summary", fun = "mean"
           )+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Month")+
  scale_x_continuous(breaks = 1:12,
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
