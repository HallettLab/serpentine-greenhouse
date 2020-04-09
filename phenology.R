######################
##Phenology clean up##
######################

##load libraries
library(tidyverse)

##datpaths
#mac
datpath <- "~/OneDrive - University Of Oregon/Sofi phenology project/"
#pc
datpath <- "C:/Users/eliza/OneDrive - University Of Oregon/Serpentine/Sofi phenology project/"

##read in data
phen_back <- read.csv(paste(datpath, "/phenology_background.csv", sep = ""))
phen_phyt <- read.csv(paste(datpath, "/phenology_phytometers.csv", sep = ""))
pots <- read.csv(paste(datpath, "/pot_treatments.csv", sep = ""))

##join trts to datasets
back_trts <- left_join(phen_back, pots)
phyt_trts <- left_join(phen_phyt, pots)

##forbs data visualization
forbs_back <- back_trts %>%
  filter(seed_sp == "LAPL" | seed_sp == "PLER")
ggplot(forbs_back) + geom_line(aes(week,flowers, color=trt_N, linetype=trt_water), size = 1) + facet_grid(seed_sp~seed_density, scale = "free")

##grasses data visualization
grass_back <- back_trts %>%
  filter(seed_sp == "BRHO" | seed_sp == "FEMI")
######################
##Phenology clean up##
######################

##load libraries
library(tidyverse)

##datpaths
#mac
datpath <- "~/OneDrive - University Of Oregon/Sofi phenology project/"
#pc
datpath <- "C:/Users/eliza/OneDrive - University Of Oregon/Serpentine/Sofi phenology project/"

##read in data
phen_back <- read.csv(paste(datpath, "/phenology_background.csv", sep = ""))
phen_phyt <- read.csv(paste(datpath, "/phenology_phytometers.csv", sep = ""))
pots <- read.csv(paste(datpath, "/pot_treatments.csv", sep = ""))

##join trts to datasets
back_trts <- left_join(phen_back, pots)
phyt_trts <- left_join(phen_phyt, pots)

##forbs data visualization
forbs_back <- back_trts %>%
  filter(seed_sp == "LAPL" | seed_sp == "PLER")
ggplot(forbs_back) + geom_point(aes(week, flowers, color=trt_N), size = 1.5) + geom_line(aes(week,flowers, color=trt_N, linetype=trt_water), size = 1) + facet_grid(seed_sp~seed_density, scale = "free") 
ggplot(forbs_back) + geom_point(aes(week, stalks_buds, color=trt_N), size = 1.5) + geom_line(aes(week,flowers, color=trt_N, linetype=trt_water), size = 1) + facet_grid(seed_sp~seed_density, scale = "free") 


##grasses data visualization
grass_back <- back_trts %>%
  filter(seed_sp == "BRHO" | seed_sp == "FEMI")
ggplot(grass_back) + geom_point(aes(week, culms, color=trt_N), size = 1.5) + geom_line(aes(week,culms, color=trt_N, linetype=trt_water), size = 1) + facet_grid(seed_sp~seed_density, scale = "free") 