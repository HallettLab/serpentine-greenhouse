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
  filter(phytometer == "PLER") %>%
  filter(!seed_sp == "none") %>%
  filter(!num_seeds =="NA")

pler_phyt2 <- seeds_phyt %>%
  filter(phytometer == "PLER") %>%
  filter(!num_seeds =="NA")

noseed_pler <- seeds_phyt %>%
  filter(phytometer == "PLER") %>%
  filter(seed_sp == "none")
  filter(!num_seeds == "NA")

p1 <- ggplot(pler_phyt) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")

p2 <- ggplot(noseed_pler) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water))

plot_grid(p1,p2)
ggsave("pler_phyt_seeds.png")

##BRHO phytometers
brho_phyt <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(!seed_sp == "none") %>%
  filter(!num_seeds == "NA")

## graph to see if viable seed changes with treatment
brho_phyt2 <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(!num_seeds == "NA") %>%
  mutate(p_viable = num_seeds/potential_seed) %>%
  filter(!p_viable == "NaN")

ggplot(brho_phyt2) + geom_boxplot(aes(trt_N, p_viable, fill = trt_water)) +
  facet_grid(~block)

ggplot(brho_phyt2) + geom_boxplot(aes(trt_N, p_viable, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")


noseed_brho <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(seed_sp == "none") %>%
  filter(!num_seeds == "NA")

b1a <- ggplot(brho_phyt) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")

b1b <- ggplot(brho_phyt) + geom_boxplot(aes(trt_N, potential_seed, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")

b2a <- ggplot(noseed_brho) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water))

b2b <- ggplot(noseed_brho) + geom_boxplot(aes(trt_N, potential_seed, fill = trt_water))

plot_grid(b1a,b2)
ggsave("brho_phyt_seeds.png")

plot_grid(b1b,b2b)
ggsave("brho_phyt_potseeds.png")

##LAPL phytometers
lapl_phyt <- seeds_phyt %>%
  filter(phytometer == "LAPL") %>%
  filter(!seed_sp == "none") %>%
  filter(!num_seeds =="NA")

lapl_phyt2 <- seeds_phyt %>%
  filter(phytometer == "LAPL") %>%
  filter(!num_seeds =="NA")

noseed_lapl <- seeds_phyt %>%
  filter(phytometer == "LAPL") %>%
  filter(seed_sp == "none") %>%
  filter(!num_seeds == "NA")

l1a <- ggplot(lapl_phyt) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")

l2a <- ggplot(noseed_lapl) + geom_boxplot(aes(trt_N, num_seeds, fill = trt_water))

l1b <- ggplot(lapl_phyt) + geom_boxplot(aes(trt_N, potential_seed, fill = trt_water)) +
  facet_grid(seed_sp~seed_density, scale = "free")

l2b <- ggplot(noseed_lapl) + geom_boxplot(aes(trt_N, potential_seed, fill = trt_water))

plot_grid(l1a,l2a)
ggsave("lapl_phyt_seeds.png")

plot_grid(l1b,l2b)
ggsave("lapl_phyt_potseeds.png")

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

plantago.labs <- c("Hi Plantago density", "Lo Plantago density", "No seed addition")
names(plantago.labs) <- c("hi", "lo", "none")

##Plots by species competition
ggplot(pler_brho) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = brho.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("pler_brho_seeds.png")

ggplot(pler_femi) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = femi.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("pler_femi_seeds.png")

ggplot(pler_lapl) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = lapl.labs)) + 
  xlab("Nitrogen trt") + ylab("# Plantago seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("pler_lapl_seeds.png")

##BRHO in competition subsets
brho_pler <- brho_phyt2 %>%
  filter(seed_sp == "PLER" | seed_sp == "none")

brho_femi <- brho_phyt2 %>%
  filter(seed_sp == "FEMI" | seed_sp == "none")

brho_lapl <- brho_phyt2 %>%
  filter(seed_sp == "LAPL" | seed_sp == "none")

##Plots by species competition
ggplot(brho_pler) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = plantago.labs)) + 
  xlab("Nitrogen trt") + ylab("# Bromus seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("brho_pler_seeds.png")

ggplot(brho_femi) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = femi.labs)) + 
  xlab("Nitrogen trt") + ylab("# Bromus seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("brho_femi_seeds.png")

ggplot(brho_lapl) + geom_boxplot(aes(trt_N, num_seeds, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = lapl.labs)) + 
  xlab("Nitrogen trt") + ylab("# Bromus seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("brho_lapl_seeds.png")

##LAPL in competition subsets
lapl_pler <- lapl_phyt2 %>%
  filter(seed_sp == "PLER" | seed_sp == "none")

lapl_femi <- lapl_phyt2 %>%
  filter(seed_sp == "FEMI" | seed_sp == "none")

lapl_brho <- lapl_phyt2 %>%
  filter(seed_sp == "BRHO" | seed_sp == "none")

##Plots by species competition
ggplot(lapl_pler) + geom_boxplot(aes(trt_N, potential_seed, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = plantago.labs)) + 
  xlab("Nitrogen trt") + ylab("# Layia seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("lapl_pler_potseeds.png")

ggplot(lapl_femi) + geom_boxplot(aes(trt_N, potential_seed, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = femi.labs)) + 
  xlab("Nitrogen trt") + ylab("# Layia seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("lapl_femi_potseeds.png")

ggplot(lapl_brho) + geom_boxplot(aes(trt_N, potential_seed, fill=trt_water)) + 
  facet_grid(~seed_density, labeller = labeller(seed_density = brho.labs)) + 
  xlab("Nitrogen trt") + ylab("# Layia seeds") + theme_bw() +
  scale_fill_discrete(name = "Water trt")
ggsave("lapl_brho_potseeds.png")