#####################
######Phenology######
#####################

##load libraries
library(tidyverse)

##datpaths
#mac Sofi
datpath <- "~/OneDrive - University Of Oregon/Sofi phenology project/"
#mac other
datpath <- "~/OneDrive - University Of Oregon/Serpentine/Sofi phenology project/"
#pc Eliza
datpath <- "C:/Users/eliza/OneDrive - University Of Oregon/Serpentine/Sofi phenology project/"

##read in data
phen_back <- read.csv(paste(datpath, "/phenology_background.csv", sep = ""))
phen_phyt <- read.csv(paste(datpath, "/phenology_phytometers.csv", sep = ""))
pots <- read.csv(paste(datpath, "/pot_treatments.csv", sep = ""))

##join trts to datasets
back_trts <- left_join(phen_back, pots) %>%
  filter(!week == 13)
phyt_trts <- left_join(phen_phyt, pots) %>%
  filter(!week == 13)

##forbs visualization
forbs_back <- back_trts %>%
  filter(seed_sp == "LAPL" | seed_sp == "PLER") 
ggplot(forbs_back) + geom_point(aes(week, flowers, color=trt_N), size = 1.5) +
  geom_line(aes(week,flowers, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free") 
ggplot(forbs_back) + geom_point(aes(week, stalks_buds, color=trt_N), size = 1.5) + 
  geom_line(aes(week,stalks_buds, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free") 

##grasses visualization
grass_back <- back_trts %>%
  filter(seed_sp == "BRHO" | seed_sp == "FEMI")
ggplot(grass_back) + geom_point(aes(week, culms, color=trt_N), size = 1.5) + 
  geom_line(aes(week,culms, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free") 

##PLER phyto visualization
pler_phyt <- phyt_trts %>%
  filter(phytometer == "PLER") 
pler_phyt[is.na(pler_phyt)] = 0
ggplot(pler_phyt) + geom_point(aes(week, flowers, color=trt_N), size = 1.5) +
  geom_line(aes(week,flowers, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free") 
ggplot(pler_phyt) + geom_point(aes(week, stalks_buds, color=trt_N), size = 1.5) +
  geom_line(aes(week,stalks_buds, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free") 

##LAPL phyto visualization
lapl_phyt <- phyt_trts %>%
  filter(phytometer == "LAPL") 
lapl_phyt[is.na(lapl_phyt)] = 0
ggplot(lapl_phyt) + geom_point(aes(week, flowers, color=trt_N), size = 1.5) +
  geom_line(aes(week,flowers, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free") 
ggplot(lapl_phyt) + geom_point(aes(week, stalks_buds, color=trt_N), size = 1.5) +
  geom_line(aes(week,stalks_buds, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free")

##BRHO phyto visualization
brho_phyt <- phyt_trts %>%
  filter(phytometer == "BRHO") 
brho_phyt[is.na(brho_phyt)] = 0
ggplot(brho_phyt) + geom_point(aes(week, culms, color=trt_N), size = 1.5) +
  geom_line(aes(week,culms, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free") 
 
##FEMI phyto visualization
femi_phyt <- phyt_trts %>%
  filter(phytometer == "FEMI") 
femi_phyt[is.na(femi_phyt)] = 0
ggplot(femi_phyt) + geom_point(aes(week, culms, color=trt_N), size = 1.5) +
  geom_line(aes(week,culms, color=trt_N, linetype=trt_water), size = 1) + 
  facet_grid(seed_sp~seed_density, scale = "free") 