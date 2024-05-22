rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)

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

###### main effect ######
df.annual.litter = 
  df.all.long |> 
  group_by(plotID, myc, tree_species_richness, block, composition) |>
  summarise(litter.prod = sum(dryweight, na.rm = T), 
            sd.litter.prod = sd(dryweight, na.rm = T))

main.eff <- ggplot(df.annual.litter, 
                   aes(x=tree_species_richness, y=litter.prod, 
                       color = myc, fill = myc))+
  geom_point(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Total litterfall"~(g/m^2)), 
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
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(a)")
main.eff


###### monthly effect ######

df.monthly.litter = df.all.long |> 
  group_by(block, plotID, sr, div, myc, month1) |> 
  summarise(litter.prod = sum(dryweight, na.rm = T))

df.monthly.litter$month1 = factor(df.monthly.litter$month1, levels = 
                                c( month.name[3:12],  month.name[1:2]), labels=month.abb)

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


fig.month <- ggplot()+
  geom_point(data = df.monthly.litter, 
             aes(x=sr, y=litter.prod, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.monthly.litter, 
              aes(x=sr, y=litter.prod, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Monthly litterfall"~(g/m^2)), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_y_continuous(limits = c(0, 160))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide = "none")+
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr'), 
            aes(x=1, y=155,
                label = sign,
                hjust = 0), 
            color = 'black') + 
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'myc'), 
            aes(x=1, y=140,
                label = sign,
                hjust = 0), 
            color = 'black') + 
  geom_text(data = M |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr:myc'), 
            aes(x=1, y=130,
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
        axis.title.x = element_text(size = 12),
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
  labs(tag = "(b)")
fig.month

# ggsave("3-plots/2-2-4-Figure-monthly-effects-sig-2024-05-07.jpeg", 
#        fig.month, 
#        height=16,
#        width=34, 
#        unit="cm", 
#        dpi=2000) 
# 
# ggsave("3-plots/2-2-4-Figure-monthly-effects-sig-2024-05-07.pdf", 
#        fig.month,
#        device = cairo_pdf,
#        height=16,
#        width=34, 
#        unit="cm", 
#        dpi=2000)


##### cumulative effect ######

df.cum.litter = df.all.long |> 
  group_by(block, plotID, div, sr, myc, month1) |> 
  summarise(litterfall = sum(dryweight, na.rm = T))

df.cum.litter$month1 = factor(df.cum.litter$month1, levels = 
                                c( month.name[3:12],  month.name[1:2]), labels=month.abb)
df.cum.litter.1 = 
  df.cum.litter |> 
  group_by(block, plotID, div, sr, myc) |> 
  arrange(month1) |> 
  mutate(cs = cumsum(litterfall))

Mc = map_df( .x = unique(df.cum.litter.1$month1),
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
Mc

Mc$sign[Mc$explanatory == 'sr' & Mc$month == 'March'] = 
  paste0("sr = ",Mc$sign[Mc$explanatory == 'sr' & Mc$month == 'March'])
Mc$sign[Mc$explanatory == 'myc' & M$month == 'March'] = 
  paste0("myc = ",Mc$sign[Mc$explanatory == 'myc' & Mc$month == 'March'])
Mc$sign[Mc$explanatory == 'sr:myc' & Mc$month == 'March'] = 
  paste0("sr:myc = ",Mc$sign[Mc$explanatory == 'sr:myc' & Mc$month == 'March'])



cum_facet<- ggplot()+
  geom_point(data = df.cum.litter.1, 
             aes(x=sr, y=cs, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.cum.litter.1, 
              aes(x=sr, y=cs, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Cumulative sum"~(g/m^2)), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_y_continuous(limits = c(0, 350))+
  #coord_cartesian(ylim = c(0, 140), xlim = c(1,4))+ 
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  geom_text(data = Mc |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr'), 
            aes(x=1, y=345,
                label = sign,
                hjust = 0), 
            color = 'black') + 
  geom_text(data = Mc |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'myc'), 
            aes(x=1, y=325,
                label = sign,
                hjust = 0), 
            color = 'black') + 
  geom_text(data = Mc |> 
              mutate(month1 = month) |> 
              filter(explanatory == 'sr:myc'), 
            aes(x=1, y=305,
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
cum_facet

# ggsave("3-plots/2-2-4-2-Figure-cumulative-sum-monthly-sig-2024-05-07.jpeg", 
#        cum_facet, 
#        height=16,
#        width=34, 
#        unit="cm", 
#        dpi=2000) 
# 
# ggsave("3-plots/2-2-4-2-Figure-cumulative-sum-monthly-sig-2024-05-07.pdf", 
#        cum_facet,
#        device = cairo_pdf,
#        height=16,
#        width=34, 
#        unit="cm", 
#        dpi=2000) 


library(gridExtra)
library(ggpubr)

ggarrange(
  main.eff,
  ggarrange(fig.month, cum_facet, ncol = 1),
  nrow = 1, ncol = 2, 
  widths = c(.3,.6)
)

plot1 <-grid.arrange(layout_matrix = rbind(c(1,2),
                                           c(1,3)),
                             grobs= list(main.eff, fig.month, cum_facet))

plot1

ggsave("3-plots/2-2-11-Figure-total-monthly-cumulative-2024-05-21.jpeg",
       plot1,
       height=28,
       width=24,
       unit="cm",
       dpi=2000)

ggsave("3-plots/2-2-11-Figure-total-monthly-cumulative-2024-05-21.pdf",
       plot1,
       device = cairo_pdf,
       height=20,
       width=34,
       unit="cm",
       dpi=2000)

