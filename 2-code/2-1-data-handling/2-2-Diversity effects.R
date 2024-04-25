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
