#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-04-25
# diversity effects 
# by Elisabeth Bönisch (elisabeth.boenisch@idiv.de)
# updated on 

#============================ Packages ===============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)

library(devtools)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
load_all("rBExIS")
bexis.options("base_url" = "https://mydivdata.idiv.de")
bexis.get.datasets()


df.all.wide.info <- read.csv("2-1-1-Full-data-wideformat-MyDiv-litter-dryweight.csv")

df.all.long = read_csv('1-data/2-1-data-handling/2-1-2-Full-data-longformat_corrected-MyDiv-litter-dryweight.csv')

df = df.all.long |> 
  group_by(block, plotID, div, myc, month1) |> 
  summarise(litterfall = sum(dryweight, na.rm = T))
df$month1 = factor(df$month1, levels = 
                     c( month.name[3:12],  month.name[1:2]))
df.1 = 
  df |> 
  group_by(block, plotID, div, myc) |> 
  arrange(month1) |> 
  mutate(cs = cumsum(litterfall))


ggplot(data = df.1, aes(x = as.numeric(month1), y = cs, color = factor(myc))) + 
  geom_smooth()
