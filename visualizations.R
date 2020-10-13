##########################
######Visualization#######
##########################

library(tidyverse)
library(minpack.lm)
library(nlstools)
library(grid)
library(gridExtra)

## Read in data
seed_dat <- read.csv(paste(datpath, "/clean_seed_dat.csv", sep = "")) 

## Data manipullation
dat <- seed_dat %>%
  group_by(block,seed_density,individual,background_comp,trt_water,trt_N) %>%
  summarize(mean_out_in = mean(out_in))

ggplot(subset(dat, individual %in% "BRHO")) + geom_boxplot(aes(trt_N,mean_out_in,fill=trt_water)) +
  facet_grid(seed_density~background_comp, scale = "free")

ggplot(subset(dat, individual %in% "PLER")) + geom_boxplot(aes(trt_N,mean_out_in,fill=trt_water)) +
  facet_grid(seed_density~background_comp, scale = "free")

ggplot(subset(dat, individual %in% "FEMI")) + geom_boxplot(aes(trt_N,mean_out_in,fill=trt_water)) +
  facet_grid(seed_density~background_comp, scale = "free")

ggplot(subset(dat, individual %in% "FEMI")) + geom_boxplot(aes(trt_N,mean_out_in,fill=trt_water)) +
  facet_grid(seed_density~background_comp, scale = "free")

dat_excnone <- dat %>%
  filter(background_comp != "none")

ggplot(subset(dat_excnone, individual %in% c("BRHO","PLER")), aes(background_comp,mean_out_in,group = interaction(individual, trt_water))) +
  geom_point(aes(color = individual, shape = trt_water),size=2.5) + geom_line(aes(color = individual)) +
  facet_grid(seed_density~trt_N) + ylab("seeds per individual") + xlab("background competition") +
  scale_shape_manual(values = c(1,16))



