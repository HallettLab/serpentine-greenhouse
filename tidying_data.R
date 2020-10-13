##############################################
##########Data clean up for analysis##########
##############################################

library(tidyverse)

## Read in data
femi_back <- read.csv(paste(datpath, "/femi_back.csv", sep = "")) 
brho_back <- read.csv(paste(datpath, "/brho_back.csv", sep = ""))
brho_phyt <- read.csv(paste(datpath, "/brho_phyt.csv", sep = ""))
back <- read.csv(paste(datpath, "/seed_biomass_background.csv", sep = "")) 
phyt <- read.csv(paste(datpath, "/seed_biomass_phytometers.csv", sep = ""))

## BRHO phytometer and background seed/biomass production
brho_back <- brho_back %>%
  rename(individual = seed_sp) %>%
  mutate(background_comp = "BRHO") %>%
  select(-X)

brho_phyt <- brho_phyt %>%
  rename(individual = phytometer, background_comp = seed_sp) %>%
  select(-X) %>%
  mutate(seeds_in = 1, seeds_out = out_in)

## FEMI phytometer and background seed/biomass production
femi_back <- femi_back %>%
  rename(individual = seed_sp) %>%
  mutate(background_comp = "FEMI") %>%
  select(-X)

femi_phyt <- phyt %>%
  filter(phytometer == "FEMI") %>%
  rename(out_in = seeds_mat) %>%
  select(-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes) %>%
  rename(individual = phytometer, background_comp = seed_sp) %>%
  mutate(seeds_in = 1, seeds_out = out_in)

## PLER phytometer and background seed/biomass production
pler_back <- back %>%
  filter(seed_sp == "PLER") %>%
  rename(seeds_out = seeds_mat) %>%
  select(-seeds_immat,-empty_glumes,-seedmass_g,-glumemass_g,-nodes,-senescence) %>%
  mutate(out_in = seeds_out/seeds_in) %>%
  rename(individual = seed_sp) %>%
  mutate(background_comp = "PLER")

pler_phyt <- phyt %>%
  filter(phytometer == "PLER") %>%
  rename(out_in = seeds_mat) %>%
  select(-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes) %>%
  rename(individual = phytometer, background_comp = seed_sp) %>%
  mutate(seeds_in = 1, seeds_out = out_in)

## LAPL phytometer and and background seed/biomass production
lapl_back <- back %>%
  filter(seed_sp == "LAPL") %>%
  mutate(seeds_out = seeds_mat + seeds_immat, out_in = seeds_out/seeds_in) %>%
  select(-seeds_mat,-seeds_immat,-empty_glumes,-seedmass_g,-glumemass_g,-nodes,-senescence) %>%
  rename(individual = seed_sp) %>%
  mutate(background_comp = "LAPL")

lapl_phyt <- phyt %>%
  filter(phytometer == "LAPL") %>%
  mutate(out_in = seeds_mat + seeds_immat) %>%
  select(-seeds_mat,-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes) %>%
  rename(individual = phytometer, background_comp = seed_sp) %>%
  mutate(seeds_in = 1, seeds_out = out_in)

## Join all species
brho_dat <- full_join(brho_back,brho_phyt)

femi_dat <- full_join(femi_back,femi_phyt)

pler_dat <- full_join(pler_back,pler_phyt)

lapl_dat <- full_join(lapl_back,lapl_phyt)

join1 <- full_join(brho_dat,femi_dat)

join2 <- full_join(pler_dat,lapl_dat)

seed_dat <- full_join(join1,join2)


