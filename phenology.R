#####################
######Phenology######
#####################

##load libraries
library(tidyverse)

#SE function
calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

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

##forbs data manipulation
forbs_back <- back_trts %>%
  filter(seed_sp == "LAPL" | seed_sp == "PLER") %>%
  group_by(week, seed_sp, seed_density, trt_water, trt_N) %>%
  summarize(avgflowers = mean(flowers), avgbuds = mean(stalks_buds), 
            seflowers = calcSE(flowers), sebuds = calcSE(stalks_buds))

##forbs data visualization
week <- 1:12
ggplot(forbs_back) + geom_point(aes(week, avgflowers, color=trt_N), size = 1.5) +
  geom_line(aes(week,avgflowers, color=trt_N, linetype=trt_water), size = 1) +
  geom_errorbar(aes(x = week, y = avgflowers, ymin = avgflowers - seflowers, ymax = avgflowers + seflowers, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)
 
ggplot(forbs_back) + geom_point(aes(week, avgbuds, color=trt_N), size = 1.5) + 
  geom_line(aes(week,avgbuds, color=trt_N, linetype=trt_water), size = 1) +
  geom_errorbar(aes(x = week, y = avgbuds, ymin = avgbuds - sebuds, ymax = avgbuds + sebuds, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)
  
##grasses data manipulation
grass_back <- back_trts %>%
  filter(seed_sp == "BRHO" | seed_sp == "FEMI")
grass_back[is.na(grass_back)] <- 0
grass_back <- grass_back %>%
  group_by(week, seed_sp, seed_density, trt_water, trt_N) %>%
  summarize(avgculms = mean(culms), seculms = calcSE(culms))

##grasses data visualization
ggplot(grass_back) + geom_point(aes(week, avgculms, color=trt_N), size = 1.5) + 
  geom_line(aes(week,avgculms, color=trt_N, linetype=trt_water), size = 1) + 
  geom_errorbar(aes(x = week, y = avgculms, ymin = avgculms - seculms, ymax = avgculms + seculms, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)

##PLER phyto data manipulation
pler_phyt <- phyt_trts %>%
  filter(phytometer == "PLER")
pler_phyt[is.na(pler_phyt)] = 0
pler_phyt <- pler_phyt %>%
  group_by(week, seed_sp, seed_density, trt_water, trt_N) %>%
  summarize(avgflowers = mean(flowers), avgbuds = mean(stalks_buds), seflowers = calcSE(flowers), sebuds = calcSE(stalks_buds))

##PLER phyto data visualization
ggplot(pler_phyt) + geom_point(aes(week, avgflowers, color=trt_N), size = 1.5) +
  geom_line(aes(week,avgflowers, color=trt_N, linetype=trt_water), size = 1) +
  geom_errorbar(aes(x = week, y = avgflowers, ymin = avgflowers - seflowers, ymax = avgflowers + seflowers, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)

ggplot(pler_phyt) + geom_point(aes(week, avgbuds, color=trt_N), size = 1.5) +
  geom_line(aes(week,avgbuds, color=trt_N, linetype=trt_water), size = 1) + 
  geom_errorbar(aes(x = week, y = avgbuds, ymin = avgbuds - sebuds, ymax = avgbuds + sebuds, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)

##LAPL phyto data manipulation
lapl_phyt <- phyt_trts %>%
  filter(phytometer == "LAPL")
lapl_phyt[is.na(lapl_phyt)] = 0
lapl_phyt <- lapl_phyt %>%
  group_by(week, seed_sp, seed_density, trt_water, trt_N) %>%
  summarize(avgflowers = mean(flowers), avgbuds = mean(stalks_buds), seflowers = calcSE(flowers), sebuds = calcSE(stalks_buds))

##LAPL phyto data visualization
ggplot(lapl_phyt) + geom_point(aes(week, avgflowers, color=trt_N), size = 1.5) +
  geom_line(aes(week, avgflowers, color=trt_N, linetype=trt_water), size = 1) + 
  geom_errorbar(aes(x = week, y = avgflowers, ymin = avgflowers - seflowers, ymax = avgflowers + seflowers, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)

ggplot(lapl_phyt) + geom_point(aes(week, avgbuds, color=trt_N), size = 1.5) +
  geom_line(aes(week, avgbuds, color=trt_N, linetype=trt_water), size = 1) +
  geom_errorbar(aes(x = week, y = avgbuds, ymin = avgbuds - sebuds, ymax = avgbuds + sebuds, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)

##BRHO phyto data manipulation
brho_phyt <- phyt_trts %>%
  filter(phytometer == "BRHO") 
brho_phyt[is.na(brho_phyt)] = 0
brho_phyt <- brho_phyt %>%
  group_by(week, seed_sp, seed_density, trt_water, trt_N) %>%
  summarize(avgculms = mean(culms), seculms = calcSE(culms))

##BRHO phyto data visualization
ggplot(brho_phyt) + geom_point(aes(week, avgculms, color=trt_N), size = 1.5) +
  geom_line(aes(week, avgculms, color=trt_N, linetype=trt_water), size = 1) + 
  geom_errorbar(aes(x = week, y = avgculms, ymin = avgculms - seculms, ymax = avgculms + seculms, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)
 
##FEMI phyto data manipulation
femi_phyt <- phyt_trts %>%
  filter(phytometer == "FEMI") 
femi_phyt[is.na(femi_phyt)] = 0
femi_phyt <- femi_phyt %>%
  group_by(week, seed_sp, seed_density, trt_water, trt_N) %>%
  summarize(avgculms = mean(culms), seculms = calcSE(culms))

##FEMI phyto data visualizaion
ggplot(femi_phyt) + geom_point(aes(week, avgculms, color=trt_N), size = 1.5) +
  geom_line(aes(week, avgculms, color=trt_N, linetype=trt_water), size = 1) +
  geom_errorbar(aes(x = week, y = avgculms, ymin = avgculms - seculms, ymax = avgculms + seculms, color = trt_N, linetype = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("week", labels = as.character(week), breaks = week)
