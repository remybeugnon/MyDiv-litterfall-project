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
  summarise(litterfall_sd = mean(litterfall, na.rm=T)/sd(litterfall, na.rm=T))|>  # mean between traps and sd between traps
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID) |> 
  summarise(litterfall_sd = mean(litterfall_sd, na.rm=T))

# use log and variability (1/stability)
hist(log(1/df.year.spat$litterfall_sd))

df.year.spat.no50 = df.year.spat[!(df.year.spat$plotID == 50),]

# 3.2) Plot: yearly spatial stability ####
spa.stab.yearly <- ggplot()+
  geom_point(data = df.year.spat , 
             aes(x=sr, y=1/litterfall_sd, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.spat, 
              aes(x=sr, y=1/litterfall_sd, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Yearly spatial variability"), 
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

#ggsave("3-plots/2-2-10-Figure-litterfalll-spatial-variability-yearly-2024-05-13.jpeg", 
       spa.stab.yearly, 
       height=16,
       width=20, 
       unit="cm", 
       dpi=2000) 

spa.stab.yearly.no50 <- ggplot()+
  geom_point(data = df.year.spat.no50 , 
             aes(x=sr, y=litterfall_sd, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.spat.no50, 
              aes(x=sr, y=litterfall_sd, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Yearly spatial variability"), 
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

ggsave("3-plots/2-2-10-Figure-litterfalll-spatial-variability-yearly-noplot50-2024-05-13.jpeg", 
       spa.stab.yearly.no50, 
       height=16,
       width=20, 
       unit="cm", 
       dpi=2000) 

#3.3) Model ####
#library(nlme)
mod.yearly.var.spat =
  lmerTest::lmer(1/litterfall_sd~ sr * myc + 
                   (1|block),
                 data = df.year.spat)

# 3.4) Check the model quality ####
#library(performance)
png("3-plots/2-2-10-Check-model-yearly-spatial-variability-m2-2024-05-23.png", 
    width=1000, height=1000)
performance::check_model(mod.yearly.var.spat)
dev.off()

# 3.5) Summary ####
summary(mod.yearly.var.spat)

# 3.6) Anova (Type III SS) ####
anova(mod.yearly.var.spat)






#### 4.1) Monthly spatial stability ####
df.month.spat = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_sd = mean(litterfall, na.rm=T)/sd(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]),labels=month.abb))

# 4.2) Plot: Monthly spatial stability ####
spa.stab.monthly <- ggplot()+
  geom_point(data = df.month.spat, 
             aes(x=sr, y=1/litterfall_sd, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.month.spat, 
              aes(x=sr, y=1/litterfall_sd, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Monthly spatial variability"), 
       x = "Tree species richness")+
  scale_y_log10(breaks=c(0.0001,0.001,0.01,0.1,1,10),labels=c(0.0001,0.001,0.01,0.1,1,10))+
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
            aes(x=1, y=20,
                label = sign,
                hjust = 0),
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'myc'),
            aes(x=1, y=15,
                label = sign,
                hjust = 0),
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'sr:myc'),
            aes(x=1, y=10,
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

#ggsave("3-plots/2-2-10-Figure-litterfalll-spatial-variability-monthly-2024-05-16.jpeg", 
       spa.stab.monthly, 
       height=16,
       width=34, 
       unit="cm", 
       dpi=2000) 

# 4.3) Model ####
# with correlation structure
#library(nlme)

hist(log(1/df.month.spat$litterfall_sd))

mod.monthly.stab.spat =
  lme(log(1/litterfall_sd) ~ month1 * sr * myc,
      random = ~1|block,
      data = df.month.spat,
      correlation=corCAR1(), na.action = na.omit)

# 4.4) Check the model quality ####
#library(performance)
png("3-plots/2-2-10-Check-model-monthly-spatial-variability-2024-05-16.png", 
    width=1000, height=1000)
performance::check_model(mod.monthly.stab.spat)
dev.off()

# 4.5) Summary ####
summary(mod.monthly.stab.spat)

# 4.6) Anova (Type III SS) ####
anova(mod.monthly.stab.spat)

# 4.7) Check months individually ####
M = map_df( .x = unique(df.month.spat$month1),
            .f = ~ {
              mod = 
                lme(log(1/litterfall_sd) ~ sr * myc, 
                    random= ~1|block,
                    data= df.month.spat |>
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





#### 5.1) Yearly temporal variability  ####
df.year.temp = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_mean = mean(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID) |> 
  summarise(litterfall_sd = mean(litterfall_mean, na.rm=T)/sd(litterfall_mean, na.rm=T))

# 5.2) Plot: Monthly spatial variability ####
temp.var.yearly <- ggplot()+
  geom_point(data = df.year.temp , 
             aes(x=sr, y=litterfall_sd, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.temp, 
              aes(x=sr, y=litterfall_sd, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Yearly temporal variability"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
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
temp.var.yearly

#ggsave("3-plots/2-2-10-Figure-litterfall-temporal-variability-yearly-2024-05-16.jpeg", 
       temp.stab.yearly, 
       height=16,
       width=20, 
       unit="cm", 
       dpi=2000) 

# 5.3) Model ####
#library(nlme)
mod.yearly.temp.var =
  lmerTest::lmer((litterfall_sd) ~ sr * myc + 
                   (1|block),
                 data = df.year.temp)

# 5.4) Check the model quality ####
#library(performance)
png("3-plots/2-2-10-Check-model-yearly-temporal-variability-m2-2024-05-23.png", 
    width=1000, height=1000)
performance::check_model(mod.yearly.temp.var)
dev.off()

# 5.5) Summary ####
summary(mod.yearly.temp.var)

# 5.6) Anova (Type III SS) ####
anova(mod.yearly.temp.var)



### end ###

library(gridExtra)
plot2 <-grid.arrange(layout_matrix = rbind(c(1,2),
                                           c(3,3)),
                     grobs= list(temp.stab.yearly, spa.stab.yearly, spa.stab.monthly))

plot2

ggsave("3-plots/2-2-10-Figure-spatial-temporal-variability-2024-05-21.jpeg",
       plot2,
       height=22,
       width=24,
       unit="cm",
       dpi=2000)

ggsave("3-plots/2-2-10-Figure-spatial-temporal-variability-2024-05-21.pdf",
       plot2,
       device = cairo_pdf,
       height=22,
       width=24,
       unit="cm",
       dpi=2000)
