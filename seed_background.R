#####################################
######Backgound seed production######
#####################################

## Load libraries
library(tidyverse)

## Load in data
seeds_back <- read.csv(paste(datpath, "/seed_background.csv", sep = "")) 
seeds_phyt <- read.csv(paste(datpath, "/seed_phytometers.csv", sep = ""))
stems_back <- read.csv(paste(datpath, "/stems_background.csv", sep = ""))

## Visualize PLER in/out
PLER_back <- seeds_back %>%
  filter(seed_sp == "PLER") %>%
  mutate(out_in = seeds_out/seeds_in)

ggplot(seeds_back) + geom_boxplot(aes(trt_N,out_in,fill=trt_water)) +
  ylab("Plantago seed produced/seed added") + facet_wrap(~seed_density)


## Data manipulation BRHO phytometers in PLER background
BRHO_phyt <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(seed_sp == "PLER") %>%
  filter(!num_seeds == "NA") %>%
  select(-biomass_mg,-comment)

PLER_stems <- stems_back %>%
  filter(seed_sp == "PLER") %>%
  mutate(recruit = stem_density/seed_added)

PLER_back <- inner_join(PLER_stems,PLER_back)

BRHO_no_comp_hi <- seeds_phyt %>%
  filter(seed_sp == "none") %>%
  filter(phytometer == "BRHO") %>%
  select(-biomass_mg,-comment,-seed_density) %>%
  rename(BRHO_mat_seeds = num_seeds, BRHO_pot_seeds = potential_seed) %>%
  group_by(seed_sp,trt_N,trt_water) %>%
  summarize(mature_BRHO = mean(BRHO_mat_seeds),
            all_BRHO = mean(BRHO_pot_seeds)) %>%
  pivot_longer(c("mature_BRHO","all_BRHO"),
               names_to = "type", values_to = "seeds") %>%
  mutate(seed_density = "hi")

BRHO_no_comp_lo <- seeds_phyt %>%
  filter(seed_sp == "none") %>%
  filter(phytometer == "BRHO") %>%
  select(-biomass_mg,-comment,-seed_density) %>%
  rename(BRHO_mat_seeds = num_seeds, BRHO_pot_seeds = potential_seed) %>%
  group_by(seed_sp,trt_N,trt_water) %>%
  summarize(mature_BRHO = mean(BRHO_mat_seeds),
            all_BRHO = mean(BRHO_pot_seeds)) %>%
  pivot_longer(c("mature_BRHO","all_BRHO"),
               names_to = "type", values_to = "seeds") %>%
  mutate(seed_density = "lo")

BRHO_no_comp <- full_join(BRHO_no_comp_hi,BRHO_no_comp_lo)

PLER_no_comp_lo <- seeds_phyt %>%
  filter(seed_sp == "none") %>%
  filter(phytometer == "PLER") %>%
  select(-biomass_mg,-comment,-seed_density) %>%
  rename(phyt_PLER = num_seeds) %>%
  group_by(seed_sp,trt_N,trt_water) %>%
  summarize(PLER = mean(phyt_PLER)) %>%
  pivot_longer(c("PLER"),
               names_to = "type", values_to = "seeds") %>%
  mutate(seed_density = "lo")

PLER_no_comp_hi <- seeds_phyt %>%
  filter(seed_sp == "none") %>%
  filter(phytometer == "PLER") %>%
  select(-biomass_mg,-comment,-seed_density) %>%
  rename(phyt_PLER = num_seeds) %>%
  group_by(seed_sp,trt_N,trt_water) %>%
  summarize(PLER = mean(phyt_PLER)) %>%
  pivot_longer(c("PLER"),
               names_to = "type", values_to = "seeds") %>%
  mutate(seed_density = "hi")

PLER_no_comp <- full_join(PLER_no_comp_hi,PLER_no_comp_lo)

seeds_PLER <- full_join(BRHO_phyt,PLER_back) %>%
  select(-seeds_in,-out_in,-seed_added,-recruit) %>%
  mutate(PLER_ind = seeds_out/stem_density) %>%
  rename(BRHO_mat_seeds = num_seeds, BRHO_pot_seeds = potential_seed,
         PLER_stems = stem_density, PLER_seeds = seeds_out) %>%
  group_by(seed_density,seed_sp,trt_water,trt_N) %>%
  summarize(all_BRHO = mean(BRHO_mat_seeds),
            mature_BRHO = mean(BRHO_pot_seeds),
            PLER = mean(PLER_ind)) %>%
  pivot_longer(c("all_BRHO","mature_BRHO","PLER"),
               names_to = "type", values_to = "seeds")

seeds_PLER <- full_join(seeds_PLER,BRHO_no_comp)

seeds_PLER <- full_join(seeds_PLER,PLER_no_comp)

mat_seeds_PLER <- seeds_PLER %>%
  filter(!type == "all_BRHO")

all_seeds_PLER <- seeds_PLER %>%
  filter(!type == "mature_BRHO")

ggplot(mat_seeds_PLER, aes(seed_sp,seeds)) +
  geom_point(aes(color = type, shape = trt_water)) +
  facet_grid(seed_density~trt_N) 

ggplot(all_seeds_PLER, aes(seed_sp,seeds)) +
  geom_point(aes(color = type, shape = trt_water)) +
  facet_grid(seed_density~trt_N) 






