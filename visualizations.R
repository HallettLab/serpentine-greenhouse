##########################
######Visualization#######
##########################

library(tidyverse)

## Read in data
seed_dat <- read.csv(paste(datpath, "/clean_seed_dat.csv", sep = ""))

## SE function
calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

## Data manipulation
dat <- seed_dat %>%
  group_by(block,seed_density,individual,background_comp,trt_water,trt_N) %>%
  summarize(out_in = out_in)

ggplot(subset(dat, individual %in% "BRHO")) + geom_boxplot(aes(trt_N,mean_out_in,fill=trt_water)) +
  facet_grid(seed_density~background_comp, scale = "free")

ggplot(subset(dat, individual %in% "PLER")) + geom_boxplot(aes(trt_N,mean_out_in,fill=trt_water)) +
  facet_grid(seed_density~background_comp, scale = "free")

ggplot(subset(dat, individual %in% "FEMI")) + geom_boxplot(aes(trt_N,mean_out_in,fill=trt_water)) +
  facet_grid(seed_density~background_comp, scale = "free")

ggplot(subset(dat, individual %in% "FEMI")) + geom_boxplot(aes(trt_N,mean_out_in,fill=trt_water)) +
  facet_grid(seed_density~background_comp, scale = "free")

dat_exc_none <- dat %>%
  group_by(seed_density,individual,background_comp,trt_water,trt_N) %>%
  summarize(mean_out_in = mean(out_in, na.rm = T),se_out_in = calcSE(out_in)) %>%
  filter(background_comp != "none")

dat_exc_none$trt_N <- factor(dat_exc_none$trt_N , levels = c("lo","int","hi"))

ggplot(dat_exc_none, aes(trt_N,mean_out_in,group = interaction(individual, trt_water))) +
  geom_point(aes(color = individual, shape = trt_water),size=2.5) + geom_line(aes(color = individual)) +
  facet_grid(seed_density~background_comp, scale = "free") + ylab("seed production per capita") + xlab("N treatments") +
  scale_shape_manual(values = c(1,16))

ggplot(dat_exc_none, aes(trt_N,mean_out_in,group = interaction(individual, trt_water))) +
  geom_point(aes(color = individual, shape = trt_water),size=2.5) + geom_line(aes(color = individual)) +
  facet_grid(seed_density~background_comp, scale = "free") + ylab("seed production per capita") + xlab("N treatments") +
  scale_shape_manual(values = c(1,16)) + 
  geom_errorbar(aes(x = trt_N, y = mean_out_in, ymin = mean_out_in - se_out_in, ymax = mean_out_in + se_out_in, 
                    color = individual), width = 0.2)



