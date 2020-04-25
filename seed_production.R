############################
######Phytometer seeds######
############################

##Load libraries
library(tidyverse)
library(cowplot)

##read in data
seeds_phyt <- read.csv(paste(datpath, "/seed_phytometers.csv", sep = ""))

##PLER phytometers
pler_phyt <- seeds_phyt %>%
  filter(phytometer == "PLER")

##BRHO phytometers
brho_phyt <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(!num_seeds == "NA")

ggplot(brho_phyt) + geom_boxplot(aes(trt_N, potential_seed, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")



##PLER in competition subsets
pler_brho <- pler_phyt %>%
  filter(seed_sp == "BRHO" | seed_sp == "none")

pler_femi <- pler_phyt %>%
  filter(seed_sp == "FEMI" | seed_sp == "none")

pler_lapl <- pler_phyt %>%
  filter(seed_sp == "LAPL" | seed_sp == "none")

##Names for facet headings
brho.labs <- c("Hi Bromus density", "Lo Bromus density", "No seed addition")
names(brho.labs) <- c("hi", "lo", "none")

femi.labs <- c("Hi Festuca density", "Lo Festuca density", "No seed addition")
names(femi.labs) <- c("hi", "lo", "none")

lapl.labs <- c("Hi Layia density", "Lo Layia density", "No seed addition")
names(lapl.labs) <- c("hi", "lo", "none")

##Plots by species competition
ggplot(pler_brho) + geom_boxplot(aes(trt_N,seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = brho.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

ggplot(pler_femi) + geom_boxplot(aes(trt_N,seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = femi.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

ggplot(pler_lapl) + geom_boxplot(aes(trt_N,seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = lapl.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")



