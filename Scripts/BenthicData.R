# Code for analysis of benthic community  #####
# Edited by Nyssa Silbiger
# Edited on 5/18/2022

#Read in required libraries
library(tidyverse)
library(here)


#Load benthic data collected by Nyssa in March 2022
Benthic.Data <-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/BenthicData/BenthicData.csv")

#Clean data and summarise per location
Benthic.Info <-Benthic.Data %>% 
  select(-Date) %>% #remove date
  group_by(CowTagID) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% #sum 4 data points per location
  rowwise() %>% 
  mutate(total= sum(c_across(DeadCoral:Sand))) #calculate total number of observations for each location

#Calculate percent benthic cover
Benthic.Cover_All <-Benthic.Info %>%
  mutate(across(DeadCoral:Sand, ~ . / total*100), .keep = c("all")) %>%  #calculate % cover for benthic categories per tile
  select(-total)


#Calculate % cover for summarised categories: coral, CCA, MA, turf, rubble, sand, sponge, deadcoral
Benthic.Cover_Categories <-Benthic.Info %>% 
  select(-total) %>% 
  mutate(TotalCoral= sum(c_across(Porites:OtherCoral))) %>% 
  mutate(TotalAlgae = sum(c_across(c(Turbinaria:OtherMacro, DeadCoral, Rubble)))) %>% 
  select(-c(Porites:OtherCoral,Turbinaria:OtherMacro, DeadCoral, Rubble)) %>% 
  rowwise() %>% 
  mutate(total= sum(c_across(CCA:TotalAlgae))) %>%  #calculate total number of observations for each location
  mutate(across(CCA:TotalAlgae, ~ . / total*100), .keep = c("all")) %>%   #calculate % cover for benthic categories per tile
  select(-total) %>%
  mutate(TotalCalc = TotalCoral+CCA) # add column for total calcifiers


