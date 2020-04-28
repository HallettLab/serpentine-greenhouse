############################
######Phytometer seeds######
############################

##Load libraries
library(tidyverse)
library(cowplot)

##read in data
seeds_phyt <- read.csv(paste(datpath, "/seed_phytometers.csv", sep = ""))
no_seed_add <- seeds_phyt %>%
  filter(seed_density == "none", seed_sp == "none")

##PLER phytometers
pler_phyt <- seeds_phyt %>%
  filter(phytometer == "PLER") %>%
  filter(!seed_sp == "none") %>%
  filter(!num_seeds =="NA")

pler_phyt2 <- seeds_phyt %>%
  filter(phytometer == "PLER") %>%
  filter(!num_seeds =="NA")

noseed_pler <- no_seed_add %>%
  filter(phytometer == "PLER") %>%
  filter(!num_seeds == "NA")

p1 <- ggplot(pler_phyt) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")

p2 <- ggplot(noseed_pler) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water))

plot_grid(p1,p2)

##BRHO phytometers
brho_phyt <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(!seed_sp == "none") %>%
  filter(!num_seeds == "NA")

brho_phyt2 <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(!num_seeds == "NA")

noseed_brho <- no_seed_add %>%
  filter(phytometer == "BRHO") %>%
  filter(!num_seeds == "NA")

b1a <- ggplot(brho_phyt) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")

b1b <- ggplot(brho_phyt) + geom_boxplot(aes(trt_N, potential_seed, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")

b2 <- ggplot(noseed_brho) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water))

plot_grid(b1a,b2)

plot_grid(b1b,b2)

##PLER in competition subsets
pler_brho <- pler_phyt2 %>%
  filter(seed_sp == "BRHO" | seed_sp == "none")

pler_femi <- pler_phyt2 %>%
  filter(seed_sp == "FEMI" | seed_sp == "none")

pler_lapl <- pler_phyt2 %>%
  filter(seed_sp == "LAPL" | seed_sp == "none")

##Names for facet headings
brho.labs <- c("Hi Bromus density", "Lo Bromus density", "No seed addition")
names(brho.labs) <- c("hi", "lo", "none")

femi.labs <- c("Hi Festuca density", "Lo Festuca density", "No seed addition")
names(femi.labs) <- c("hi", "lo", "none")

lapl.labs <- c("Hi Layia density", "Lo Layia density", "No seed addition")
names(lapl.labs) <- c("hi", "lo", "none")

##Plots by species competition
ggplot(pler_brho) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = brho.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

ggplot(pler_femi) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = femi.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

ggplot(pler_lapl) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = lapl.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

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
ggplot(pler_brho) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = brho.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

ggplot(pler_femi) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = femi.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

ggplot(pler_lapl) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = lapl.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

##BRHO in competition subsets
brho_pler <- brho_phyt2 %>%
  filter(seed_sp == "PLER" | seed_sp == "none")

brho_femi <- brho_phyt2 %>%
  filter(seed_sp == "FEMI" | seed_sp == "none")

brho_lapl <- brho_phyt2 %>%
  filter(seed_sp == "LAPL" | seed_sp == "none")

##Names for facet headings
pler.labs <- c("Hi Bromus density", "Lo Bromus density", "No seed addition")
names(brho.labs) <- c("hi", "lo", "none")

##Plots by species competition
ggplot(brho_pler) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = brho.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

ggplot(brho_femi) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = femi.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

ggplot(brho_lapl) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = lapl.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")

