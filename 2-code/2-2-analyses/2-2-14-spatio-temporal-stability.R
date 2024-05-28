#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-05-11
# temporal stability and spatial stability of litter fall
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

# 1) NO! average across traps ####
df.all.wide.trap <- df.all.wide.info 

df.all.wide.trap$div<-as.factor(df.all.wide.trap$tree_species_richness)
df.all.wide.trap$blk<-as.factor(df.all.wide.trap$block)
df.all.wide.trap$myc<-as.factor(df.all.wide.trap$mycorrhizal_type)
df.all.wide.trap$sr<-df.all.wide.trap$tree_species_richness
df.all.wide.trap$sr_myc<-paste(df.all.wide.trap$sr,df.all.wide.trap$myc,sep="_")#interaction terms
df.all.wide.trap$myc <- recode_factor(df.all.wide.trap$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")
df.all.wide.trap$block <- as.factor(df.all.wide.trap$block)


# 2) long format ####
df.all.long <- df.all.wide.trap %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="litterfall")

#### 3.1) Yearly spatial stability ####
df.year.spat = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_stability = mean(litterfall, na.rm=T)/sd(litterfall, na.rm=T))|>  # mean between traps and sd between traps
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID) |> 
  summarise(litterfall_stability = mean(litterfall_stability, na.rm=T))

# use log for variability (1/stability)
# hist(log(1/df.year.spat$litterfall_stability))

df.year.spat.no50 = df.year.spat[!(df.year.spat$plotID == 50),]

# 3.2) Plot: yearly spatial stability ####
spa.stab.yearly <- ggplot()+
  geom_point(data = df.year.spat , 
             aes(x=sr, y=litterfall_stability, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.spat, 
              aes(x=sr, y=litterfall_stability, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Yearly spatial stability"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_y_log10() +
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
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
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(b)")
spa.stab.yearly

#ggsave("3-plots/2-2-14-Figure-litterfalll-spatial-stability-yearly-2024-05-24.jpeg", 
       spa.stab.yearly, 
       height=16,
       width=20, 
       unit="cm", 
       dpi=2000) 

spa.stab.yearly.no50 <- ggplot()+
  geom_point(data = df.year.spat.no50 , 
             aes(x=sr, y=litterfall_stability, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.spat.no50, 
              aes(x=sr, y=litterfall_stability, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Yearly spatial stability"), 
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

spa.stab.yearly.no50

#ggsave("3-plots/2-2-14-Figure-litterfall-spatial-stability-yearly-noplot50-2024-05-24.jpeg", 
       spa.stab.yearly.no50, 
       height=16,
       width=20, 
       unit="cm", 
       dpi=2000) 

#3.3) Model ####
#library(nlme)
mod.yearly.spa.stab =
  lmerTest::lmer(log(litterfall_stability) ~ log2(sr) * myc + 
                   (1|block),
                 data = df.year.spat)

hist(log(df.year.spat$litterfall_stability))

# 3.4) Check the model quality ####
#library(performance)
png("3-plots/2-2-14-Check-model-yearly-spatial-stability-m2-2024-05-24.png", 
    width=1000, height=1000)
performance::check_model(mod.yearly.spa.stab)
dev.off()

# 3.5) Summary ####
summary(mod.yearly.spa.stab)

# 3.6) Anova (Type III SS) ####
anova(mod.yearly.spa.stab)






#### 4.1) Monthly spatial stability ####
df.month.spat.stab = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_stability = mean(litterfall, na.rm=T)/sd(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]),labels=month.abb))

hist(log(df.month.spat.stab$litterfall_stability))

# 4.2) Plot: Monthly spatial stability ####

# 4.7) Check months individually ####
M = map_df( .x = unique(df.month.spat.stab$month1),
            .f = ~ {
              mod = 
                lme(log(litterfall_stability) ~ log2(sr) * myc, 
                    random= ~1|block,
                    data= df.month.spat.stab |>
                      filter(month1 == .x), 
                    na.action = na.omit) |> 
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

spa.stab.monthly <- ggplot()+
  geom_point(data = df.month.spat.stab, 
             aes(x=sr, y=litterfall_stability, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.month.spat.stab, 
              aes(x=sr, y=litterfall_stability, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Monthly spatial stability"), 
       x = "Tree species richness")+
  #coord_trans(y="log2")+
  # scale_y_log(
  #   #breaks=c(0.1,1,10,100,1000),labels=c(0.1,1,10,100,1000),
  #   )+
  scale_y_continuous(trans='log',
                     breaks=c(0.1,1,10,100,500), labels=c(0.1,1,10,100,500))+
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
            aes(x=1, y=110,
                label = sign,
                hjust = 0),
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'myc'),
            aes(x=1, y=90,
                label = sign,
                hjust = 0),
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'sr:myc'),
            aes(x=1, y=70,
                label = sign,
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
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(c)")

spa.stab.monthly

#ggsave("3-plots/2-2-14-Figure-litterfalll-spatial-stability-monthly-2024-05-24.jpeg", 
       spa.stab.monthly, 
       height=16,
       width=34, 
       unit="cm", 
       dpi=2000) 

# 4.3) Model ####
# with correlation structure
#library(nlme)

hist(log(df.month.spat.stab$litterfall_stability))

mod.monthly.stab.spat =
  lme(log(litterfall_stability) ~ month1 * log2(sr) * myc,
      random = ~1|block,
      data = df.month.spat.stab,
      correlation=corCAR1(), na.action = na.omit)

# 4.4) Check the model quality ####
#library(performance)
png("3-plots/2-2-14-Check-model-monthly-spatial-stability-2024-05-16.png", 
    width=1000, height=1000)
performance::check_model(mod.monthly.stab.spat)
dev.off()

# 4.5) Summary ####
summary(mod.monthly.stab.spat)

# 4.6) Anova (Type III SS) ####
anova(mod.monthly.stab.spat)





#### 5.1) Yearly temporal stability  ####
df.year.temp.stab = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_mean = mean(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID) |> 
  summarise(litterfall_stability = mean(litterfall_mean, na.rm=T)/sd(litterfall_mean, na.rm=T))

# 5.2) Plot: Monthly spatial stability ####
temp.stab.yearly <- ggplot()+
  geom_point(data = df.year.temp.stab , 
             aes(x=sr, y=litterfall_stability, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.temp.stab, 
              aes(x=sr, y=litterfall_stability, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Yearly temporal stability"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  #scale_y_continuous(trans='log')+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
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
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(a)")
temp.stab.yearly

#ggsave("3-plots/2-2-10-Figure-litterfall-temporal-stability-yearly-2024-05-16.jpeg", 
       temp.stab.yearly, 
       height=16,
       width=20, 
       unit="cm", 
       dpi=2000) 

# 5.3) Model ####
#library(nlme)
mod.yearly.temp.stab =
  lmerTest::lmer(litterfall_stability ~ log2(sr) * myc + 
                   (1|block),
                 data = df.year.temp.stab)

# 5.4) Check the model quality ####
#library(performance)
png("3-plots/2-2-10-Check-model-yearly-temporal-stability-m2-2024-05-23.png", 
    width=1000, height=1000)
performance::check_model(mod.yearly.temp.stab)
dev.off()

# 5.5) Summary ####
summary(mod.yearly.temp.stab)

# 5.6) Anova (Type III SS) ####
anova(mod.yearly.temp.stab)



### end ###

library(gridExtra)
plot2 <-grid.arrange(layout_matrix = rbind(c(1,2),
                                           c(3,3)),
                     grobs= list(temp.stab.yearly, spa.stab.yearly, spa.stab.monthly))

plot2

ggsave("3-plots/2-2-14-Figure-spatial-temporal-stability-m2-2024-05-27.jpeg",
       plot2,
       height=22,
       width=24,
       unit="cm",
       dpi=2000)

ggsave("3-plots/2-2-14-Figure-spatial-temporal-stability-m2-2024-05-27.pdf",
       plot2,
       device = cairo_pdf,
       height=22,
       width=24,
       unit="cm",
       dpi=2000)
