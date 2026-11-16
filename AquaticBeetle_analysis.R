library(readxl)
library(readr)
library(tidyverse)
library(bipartite)
library(patchwork)
library(forcats)
library(igraph)
library(DHARMa)

rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

set.seed(123)

#1.Find the distribution of individuals across time

df_24 <- read.csv("Foodweb_samples_2024.csv", sep = ";") |> 
  filter(Group == "Aquatic_predators")
df_25 <- read.csv("FoodwebSamples_2025.csv", sep = ";")|> 
  filter(Group == "Aquatic_predators")

df_tot <- rbind(df_24, df_25)

#Trun dates into numeric
#Clean names to clearly show adult or larvae
#Plot proportion of community against time

#2. 
#Check if diet differs between larval stages and between larvae and adults
#Get diet data
#Will require re-cleaning of specifically the C.dolabratus data
#We can work directly with the edgelists to figure this out
#some type of chi-sqr test for the distribution of FOO or meanRRA
