#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-30
# monthly effects 
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


# 3) monthly effect - sum per plot and month ####
df.monthly.litter = df.all.long |> 
  group_by(block, plotID, sr, div, myc, month1) |> 
  summarise(litter.prod = sum(dryweight, na.rm = T))

df.monthly.litter$month1 = factor(df.monthly.litter$month1, levels = 
                     c( month.name[3:12],  month.name[1:2]))


# 4) plot monthly effects ####
fig.month <- ggplot()+
  geom_point(data = df.monthly.litter, 
             aes(x=sr, y=litter.prod, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.monthly.litter, 
              aes(x=sr, y=litter.prod, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Leaf litter dryweight"~(g/m2)), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_y_continuous(limits = c(0, 140))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr'), 
            aes(x=1, y=135,
                label = sign,
                fontface = "bold",
                hjust = 0), 
            color = 'black') + 
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'myc'), 
            aes(x=1, y=130,
                label = sign,
                fontface = "bold",
                hjust = 0), 
            color = 'black') + 
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr:myc'), 
            aes(x=1, y=125,
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
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))

ggsave("3-plots/2-2-4-Figure-monthly-effects-sig-2024-05-07.jpeg", 
       fig.month, 
       height=20,
       width=28, 
       unit="cm", 
       dpi=2000) 

ggsave("3-plots/2-2-4-Figure-monthly-effects-sig-2024-05-07.pdf", 
       fig.month,
       device = cairo_pdf,
       height=20,
       width=28, 
       unit="cm", 
       dpi=2000) 

# 4) Model ####

# with correlation structure
#library(nlme)
mod.monthly.litterfall =
  lme(litter.prod ~ month1 * sr * myc,
      random = ~1|block,
                 data = df.monthly.litter,
                 correlation=corCAR1())

# 5) Check the model quality ####
#library(performance)
png("3-plots/2-2-6-Check-model-monthly-effects-2024-05-07.png", 
    width=1000, height=1000)
performance::check_model(mod.monthly.litterfall)
dev.off()

# 6) Summary ####
summary(mod.monthly.litterfall)

# 7) Anova (Type III SS) ####
anova(mod.monthly.litterfall)

# 8) Check months individually
M = map_df( .x = unique(df.monthly.litter$month1),
     .f = ~ {
       mod = 
         lme(litter.prod ~ sr * myc, 
             random= ~1|block,
              data= df.monthly.litter |>
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