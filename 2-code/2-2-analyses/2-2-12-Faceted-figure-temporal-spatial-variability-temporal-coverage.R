#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-05-22
# temporal and spatial variability + temporal coverage of litter fall
# by Elisabeth Bönisch (elisabeth.boenisch@idiv.de)

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

# NO! average across traps ####
df.all.wide.trap <- df.all.wide.info 

df.all.wide.trap$div<-as.factor(df.all.wide.trap$tree_species_richness)
df.all.wide.trap$blk<-as.factor(df.all.wide.trap$block)
df.all.wide.trap$myc<-as.factor(df.all.wide.trap$mycorrhizal_type)
df.all.wide.trap$sr<-df.all.wide.trap$tree_species_richness
df.all.wide.trap$sr_myc<-paste(df.all.wide.trap$sr,df.all.wide.trap$myc,sep="_")#interaction terms
df.all.wide.trap$myc <- recode_factor(df.all.wide.trap$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")
df.all.wide.trap$block <- as.factor(df.all.wide.trap$block)

df.all.long <- df.all.wide.trap %>% 
  pivot_longer(cols=c(8:17),
               names_to="species",
               values_to="litterfall")

#### 1 Yearly spatial stability ####
df.year.spat = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_sd = mean(litterfall, na.rm=T)/sd(litterfall, na.rm=T))|>  # mean between traps and sd between traps
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID) |> 
  summarise(litterfall_sd = mean(litterfall_sd, na.rm=T))

# use log and variability (1/stability)
hist(log(1/df.year.spat$litterfall_sd))

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
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_text(color="black", size = 12),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
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
  labs(tag = "(c)")
spa.stab.yearly

#### 2 Monthly spatial stability ####
df.month.spat = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_sd = mean(litterfall, na.rm=T)/sd(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]),labels=month.abb))

df.month.spat$month2 <- dplyr::recode_factor(df.month.spat$month1, 
                                                 "March"="Mar", 
                                                 "April"="Apr",
                                                 "May"="May",
                                                 "June"="Jun",
                                                 "July"="Jul",
                                                 "August"="Aug",
                                                 "September"="Sep",
                                                 "October"="Oct",
                                                 "November"="Nov",
                                                 "December"="Dec",
                                                 "January"="Jan",
                                                 "February"="Feb")

df.month.spat$month1 = factor(df.month.spat$month2, 
                                  levels = c(month.abb[3:12], month.abb[1:2]))


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
            size=3, 
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'myc'),
            aes(x=1, y=15,
                label = sign,
                hjust = 0),
            size=3, 
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'sr:myc'),
            aes(x=1, y=10,
                label = sign,
                hjust = 0),
            size=3, 
            color = 'black') +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_text(color="black", size = 12),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
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
  labs(tag = "(d)")

spa.stab.monthly


#### 3 Yearly temporal stability  ####
df.year.temp = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(litterfall_mean = mean(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID) |> 
  summarise(litterfall_sd = mean(litterfall_mean, na.rm=T)/sd(litterfall_mean, na.rm=T))

temp.stab.yearly <- ggplot()+
  geom_point(data = df.year.temp , 
             aes(x=sr, y=log(1/litterfall_sd), 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.temp, 
              aes(x=sr, y=log(1/litterfall_sd), 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Yearly temporal variability"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_y_log10()+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_text(color="black", size = 12),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
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


#### 5 temporal coverage ####

# average across traps ####
df.all.wide.mean <- df.all.wide.info %>%
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

# long format ####
df.all.long <- df.all.wide.mean %>% 
  pivot_longer(cols=c(13:22),
               names_to="species",
               values_to="dryweight")

# sum across species per plot
df.litter.sum = df.all.long |> 
  group_by(block, plotID, div, sr, myc, month1) |> 
  summarise(litterfall = sum(dryweight, na.rm = T),
            litterfall_mean = mean(dryweight, na.rm =T),
            litterfall_sd = sd(dryweight, na.rm = T))

# monthly coverage - sum per plot, if litter dryweight larger then 0 ####
df.litter.cover = df.litter.sum |>
  group_by(block, plotID, div, sr, myc) |>
  filter(litterfall >0) |>
  summarise(number_month_litterfall = n())

ggplot(df.litter.cover, aes(x=number_month_litterfall, y=div, color=myc, fill=myc))+
  geom_bar(position = "dodge", stat = "summary", aes(fill=myc, color=myc), alpha=0.5)+
  labs(y=bquote("Tree species richness"), 
       x = "Number of months of litterfall")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  scale_x_continuous(breaks= c(0:12))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_text(color="black", size = 12),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
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

temp.cov <- ggplot(df.litter.cover, aes(y=number_month_litterfall, x=sr, color=myc, fill=myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Number of months of litterfall"), 
       x = "Tree species richness")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  scale_y_continuous(breaks= c(0:12))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_text(color="black", size = 12),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
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

temp.cov


library(gridExtra)
plot2 <-grid.arrange(layout_matrix = rbind(c(1,2,3), 
                                           c(4,4,4)),
                     grobs= list(temp.stab.yearly, temp.cov,spa.stab.yearly, spa.stab.monthly))

plot2


ggsave("3-plots/2-2-12-Figure-spatial-temporal-variability-m2-2024-05-22.jpeg",
       plot2,
       height=20,
       width=32,
       unit="cm",
       dpi=2000)

ggsave("3-plots/2-2-12-Figure-spatial-temporal-variability-m2-2024-05-22.pdf",
       plot2,
       device = cairo_pdf,
       height=20,
       width=32,
       unit="cm",
       dpi=2000)
