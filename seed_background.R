#####################################
######Backgound seed production######
#####################################

## Load libraries
library(tidyverse)

## Load in data
seeds_back <- read.csv(paste(datpath, "/seed_background.csv", sep = "")) %>%
  filter(seed_sp == "PLER") %>%
  mutate(out_in = seeds_out/seeds_in)

## Visualize PLER in/out
ggplot(seeds_back) + geom_boxplot(aes(trt_N,out_in,fill=trt_water)) +
  ylab("Plantago seed produced/seed added")



