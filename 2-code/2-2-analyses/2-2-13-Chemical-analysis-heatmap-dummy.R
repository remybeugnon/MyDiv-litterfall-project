#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-05-23
# chemical composition
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

# 1) average across traps ####
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


# 2) long format ####
df.all.long <- df.all.wide.mean %>% 
  pivot_longer(cols=c(13:22),
               names_to="species",
               values_to="dryweight")


# monthly effect - sum per plot and month ####
df.monthly.litter = df.all.long |> 
  group_by(block, plotID, sr, div, myc, month1) |> 
  summarise(monthly_litterfall = sum(dryweight, na.rm = T))

df.monthly.litter$month2 <- dplyr::recode_factor(df.monthly.litter$month1, 
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

df.monthly.litter$month1 = factor(df.monthly.litter$month2, 
                                  levels = c(month.abb[3:12], month.abb[1:2]))

# cumulative sum - sum across species ####
df.cum.litter = df.all.long |> 
  group_by(block, plotID, div, sr, myc, month1) |> 
  summarise(litterfall = sum(dryweight, na.rm = T))

df.cum.litter.1 = 
  df.cum.litter |> 
  group_by(block, plotID, div, sr, myc) |> 
  arrange(month1) |> 
  mutate(cumulative_sum = cumsum(litterfall))

df.cum.litter.1$month2 <- dplyr::recode_factor(df.cum.litter.1$month1, 
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

df.cum.litter.1$month1 = factor(df.cum.litter.1$month2, 
                                  levels = c(month.abb[3:12], month.abb[1:2]))

# monthly spatial stability - sum across species ####
# NO! average across traps 
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

df.month.spat = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1) |> 
  summarise(spatial_variability = mean(litterfall, na.rm=T)/sd(litterfall, na.rm=T))|>
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


# join monthly litterfall, cumulative sum & spatial variability ####

df.m.cs <- full_join(df.monthly.litter, df.cum.litter.1, by=c("block", "plotID", "div", "sr", "myc", "month1", "month2"))

df.m.cs$blk<-as.factor(df.m.cs$block)
df.m.cs$block<-as.factor(df.m.cs$block)
df.m.cs$sr<-as.numeric(df.m.cs$sr)
df.m.cs$div<-as.factor(df.m.cs$sr)
df.m.cs$myc <- recode_factor(df.m.cs$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

df.m.cs.var <- full_join(df.m.cs, df.month.spat, by=c("block", "plotID", "div", "sr", "myc", "month1", "month2"))

# 
# df = 
#   tibble(
#     C = runif(120),
#     N = runif(120),
#     month = factor(month.name) |> rep(10)
#   ) |>
#   pivot_longer(cols = 1:2)
install.packages("RColorBrewer")
library("RcolorBrewer")

colfunc <- colorRampPalette(c("#89216B","#E94057","#F27121"))
colfunc(10)

heat <- 
  ggplot() + 
  geom_raster(data = df.m.cs.var, 
              aes(x = month1, y = "C", fill = cumulative_sum))+
  geom_raster(data = df.m.cs.var, 
              aes(x = month1, y = "N", fill = spatial_variability))+
  geom_raster(data = df.m.cs.var, 
              aes(x = month1, y = "P", fill = monthly_litterfall))+
  labs(y="Leaf litter elemental concentration")+
  scale_fill_gradientn(name = "Concentration",
                      colours = colfunc(10))+ 
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
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

ggsave("3-plots/2-2-13-Figure-elemental-concentration-2024-05-23.jpeg",
       heat,
       height=18,
       width=36,
       unit="cm",
       dpi=2000)  
  
ggsave("3-plots/2-2-13-Figure-elemental-concentration-2024-05-23.pdf",
       heat,
       device = cairo_pdf,
       height=18,
       width=36,
       unit="cm",
       dpi=2000)
