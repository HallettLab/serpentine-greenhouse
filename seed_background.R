#####################################
######Backgound seed production######
#####################################

## Load libraries
library(tidyverse)

## Load in data
seeds_back <- read.csv(paste(datpath, "/seed_background.csv", sep = "")) 
seeds_phyt <- read.csv(paste(datpath, "/seed_phytometers.csv", sep = ""))
stems_back <- read.csv(paste(datpath, "/stems_background.csv", sep = ""))

## Visualize BRHP in/out
BRHO_back <- seeds_back %>%
  filter(seed_sp == "BRHO") %>%
  mutate(out_in_mat = seeds_mat/seeds_in) %>%
  mutate(seeds_total = seeds_mat + seeds_immat) %>%
  mutate(out_in_total = seeds_total/seeds_in)

ggplot(BRHO_back) + geom_boxplot(aes(trt_N,out_in_mat,fill=trt_water)) +
  ylab("Bromus seed produced/seed added") + facet_wrap(~seed_density,labeller = labeller(seed_density = density.labs)) +
  scale_fill_manual(values=c("azure3","azure4"))

ggplot(BRHO_back) + geom_boxplot(aes(trt_N,out_in_total,fill=trt_water)) +
  ylab("Bromus seed produced/seed added") + facet_wrap(~seed_density,labeller = labeller(seed_density = density.labs)) +
  scale_fill_manual(values=c("azure3","azure4"))

ggplot(BRHO_back) + geom_boxplot(aes(trt_N,empty_glumes,fill=trt_water)) +
  ylab("Bromus empty glumes") + facet_wrap(~seed_density,labeller = labeller(seed_density = density.labs)) +
  scale_fill_manual(values=c("azure3","azure4"))

BRHO_back$seed_density <- factor(BRHO_back$seed_density , levels = c("lo","hi"))
BRHO_back$trt_N <- factor(BRHO_back$trt_N , levels = c("lo","int","hi"))
BRHO_back$trt_water <- factor(BRHO_back$trt_water , levels = c("lo","hi"))

## Visualize PLER in/out
PLER_back <- seeds_back %>%
  filter(seed_sp == "PLER") %>%
  mutate(out_in = seeds_out/seeds_in)

ggplot(PLER_back) + geom_boxplot(aes(trt_N,out_in,fill=trt_water)) +
  ylab("Plantago seed produced/seed added") + facet_wrap(~seed_density,,labeller = labeller(seed_density = density.labs)) +
  scale_fill_manual(values=c("azure3","azure4"))

PLER_back$seed_density <- factor(PLER_back$seed_density , levels = c("lo","hi"))
PLER_back$trt_N <- factor(PLER_back$trt_N , levels = c("lo","int","hi"))
PLER_back$trt_water <- factor(PLER_back$trt_water , levels = c("lo","hi"))


## Visualize recruitment of background
stems_back <- stems_back %>%
  filter(!stem_density == "NA") %>%
  filter(seed_sp == "PLER" | seed_sp == "BRHO") %>%
  mutate(recruitment = stem_density/seed_added)

ggplot(stems_back) + geom_boxplot(aes(trt_N,recruitment,fill=trt_water)) +
  facet_grid(seed_sp~seed_density, scale="free", labeller = labeller(seed_density = density.labs)) +
  scale_fill_manual(values=c("azure3","azure4"))
stems_back$seed_density <- factor(stems_back$seed_density , levels = c("lo","hi"))
stems_back$trt_N <- factor(stems_back$trt_N , levels = c("lo","int","hi"))
stems_back$trt_water <- factor(stems_back$trt_water , levels = c("lo","hi"))


density.labs <- c("low seed density", "high seed density")
names(density.labs) <- c("lo", "hi")







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


## PLER seeds out/stem density
seeds_PLER_stem <- full_join(BRHO_phyt,PLER_back) %>%
  select(-seed_added, -recruit, -seeds_in, -out_in) %>%
  mutate(PLER_ind = seeds_out/stem_density) %>%
  rename(BRHO_mat_seeds = num_seeds, BRHO_pot_seeds = potential_seed,
         PLER_stems = stem_density) %>%
  group_by(seed_density,seed_sp,trt_water,trt_N) %>%
  summarize(all_BRHO = mean(BRHO_pot_seeds),
            mature_BRHO = mean(BRHO_mat_seeds),
            PLER = mean(PLER_ind)) %>%
  pivot_longer(c("all_BRHO","mature_BRHO","PLER"),
               names_to = "type", values_to = "seeds")

seeds_PLER_stem <- full_join(seeds_PLER_stem,BRHO_no_comp)

seeds_PLER_stem <- full_join(seeds_PLER_stem,PLER_no_comp)

mat_seeds_PLERstem <- seeds_PLER_stem %>%
  filter(!type == "all_BRHO")

all_seeds_PLERstem <- seeds_PLER_stem %>%
  filter(!type == "mature_BRHO")

ggplot(mat_seeds_PLERstem, aes(seed_sp,seeds,group = interaction(type, trt_water))) +
  geom_point(aes(color = type, shape = trt_water)) + geom_line(aes(color = type)) +
  facet_grid(seed_density~trt_N, labeller = labeller(seed_density = density.labs, trt_N = n.labs)) + ylab("seeds per individual") + xlab("background competition")

ggplot(all_seeds_PLERstem, aes(seed_sp,seeds, group = interaction(type, trt_water))) +
  geom_point(aes(color = type, shape = trt_water)) + geom_line(aes(color = type)) + 
  facet_grid(seed_density~trt_N,labeller = labeller(seed_density = density.labs, trt_N = n.labs)) + ylab("seeds per individual") + xlab("background competition")

n.labs <- c("high N", "int N", "low N")
names(n.labs) <- c("hi","int", "lo")

## PLER seeds in/out
seeds_PLER_outin <- full_join(BRHO_phyt,PLER_back) %>%
  select(-seed_added, -recruit, -stem_density, -seeds_out,-seeds_in) %>%
  rename(BRHO_mat_seeds = num_seeds, BRHO_pot_seeds = potential_seed,
         PLER_outin = out_in) %>%
  group_by(seed_density,seed_sp,trt_water,trt_N) %>%
  summarize(all_BRHO = mean(BRHO_pot_seeds),
            mature_BRHO = mean(BRHO_mat_seeds),
            PLER = mean(PLER_outin)) %>%
  pivot_longer(c("all_BRHO","mature_BRHO","PLER"),
               names_to = "type", values_to = "seeds")

seeds_PLER_outin <- full_join(seeds_PLER_outin,BRHO_no_comp)

seeds_PLER_outin <- full_join(seeds_PLER_outin,PLER_no_comp)

mat_seeds_PLERoutin <- seeds_PLER_outin %>%
  filter(!type == "all_BRHO")

all_seeds_PLERoutin <- seeds_PLER_outin %>%
  filter(!type == "mature_BRHO")

ggplot(mat_seeds_PLERoutin, aes(seed_sp,seeds,group = interaction(type, trt_water))) +
  geom_point(aes(color = type, shape = trt_water),size=2.5) + geom_line(aes(color = type)) +
  facet_grid(seed_density~trt_N, labeller = labeller(seed_density = density.labs, trt_N = n.labs)) + ylab("seeds per individual") + xlab("background competition") +
  scale_shape_manual(values = c(1,16))

ggplot(all_seeds_PLERoutin, aes(seed_sp,seeds, group = interaction(type, trt_water))) +
  geom_point(aes(color = type, shape = trt_water),size=2.5) + geom_line(aes(color = type)) + 
  facet_grid(seed_density~trt_N,labeller = labeller(seed_density = density.labs, trt_N = n.labs)) + ylab("seeds per individual") + xlab("background competition") +
  scale_shape_manual(values = c(1,16))

all_seeds_PLERoutin$seed_density <- factor(all_seeds_PLERoutin$seed_density , levels = c("lo","hi"))
all_seeds_PLERoutin$trt_N <- factor(all_seeds_PLERoutin$trt_N , levels = c("lo","int","hi"))
all_seeds_PLERoutin$trt_water <- factor(all_seeds_PLERoutin$trt_water , levels = c("lo","hi"))

mat_seeds_PLERoutin$seed_density <- factor(mat_seeds_PLERoutin$seed_density , levels = c("lo","hi"))
mat_seeds_PLERoutin$trt_N <- factor(mat_seeds_PLERoutin$trt_N , levels = c("lo","int","hi"))
mat_seeds_PLERoutin$trt_water <- factor(mat_seeds_PLERoutin$trt_water , levels = c("lo","hi"))






