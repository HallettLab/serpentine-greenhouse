##############################################
##########Data clean up for analysis##########
##############################################

library(tidyverse)

## Read in data
femi_back <- read.csv(paste(datpath, "/femi_back.csv", sep = "")) 
brho_back <- read.csv(paste(datpath, "/brho_back.csv", sep = ""))
brho_phyt <- read.csv(paste(datpath, "/brho_phyt.csv", sep = ""))
brho_back2 <- read.csv(paste(datpath, "/brho_back2.csv", sep = ""))
brho_phyt2 <- read.csv(paste(datpath, "/brho_phyt2.csv", sep = ""))
back <- read.csv(paste(datpath, "/seed_biomass_background.csv", sep = "")) 
phyt <- read.csv(paste(datpath, "/seed_biomass_phytometers.csv", sep = ""))
dens <- read.csv(paste(datpath, "/seed_density.csv", sep = ""))

## BRHO phytometer and background seed/biomass production
brho_back <- brho_back %>%
  rename(species = seed_sp, BRHO_seeds_in = seeds_in) %>%
  mutate(background = "BRHO") %>%
  select(-X, -biomass_g, -out_in) %>%
  mutate(PLER_seeds_in = 0, LAPL_seeds_in = 0, FEMI_seeds_in = 0, waterN_treatment = paste(trt_water, trt_N, sep = "_"))

brho_back2 <- brho_back2 %>%
  select(-X, -biomass_g, -seeds_out) %>%
  rename(species = seed_sp, seeds_out=seeds_out_new, BRHO_seeds_in = seeds_in) %>%
  mutate(background = "BRHO") %>%
  mutate(PLER_seeds_in = 0, LAPL_seeds_in = 0, FEMI_seeds_in = 0, waterN_treatment = paste(trt_water, trt_N, sep = "_"))

brho_phyt <- brho_phyt %>%
  rename(species = phytometer, background = seed_sp, seeds_out = out_in) %>%
  select(-X,-biomass_g) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_")) %>%
  inner_join(dens)

brho_phyt2 <- brho_phyt2 %>%
  rename(species = phytometer, background = seed_sp) %>%
  select(-X,-seeds_in) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_")) %>%
  inner_join(dens)

## FEMI phytometer and background seed/biomass production
femi_back <- femi_back %>%
  rename(species = seed_sp, FEMI_seeds_in = seeds_in) %>%
  mutate(background = "FEMI", PLER_seeds_in = 0, BRHO_seeds_in = 0, LAPL_seeds_in = 0, waterN_treatment = paste(trt_water, trt_N, sep = "_")) %>%
  select(-X, -out_in, -biomass_g)

femi_phyt <- phyt %>%
  filter(phytometer == "FEMI") %>%
  rename(seeds_out = seeds_mat) %>%
  select(-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes,-biomass_g) %>%
  rename(species = phytometer, background = seed_sp) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_")) %>%
  inner_join(dens)

## PLER phytometer and background seed/biomass production
pler_back <- back %>%
  filter(seed_sp == "PLER") %>%
  rename(seeds_out = seeds_mat) %>%
  select(-seeds_immat,-empty_glumes,-seedmass_g,-glumemass_g,-nodes,-senescence,-biomass_g) %>%
  mutate(background = "PLER", BRHO_seeds_in = 0, LAPL_seeds_in = 0, FEMI_seeds_in = 0, waterN_treatment = paste(trt_water, trt_N, sep = "_")) %>%
  rename(species = seed_sp, PLER_seeds_in = seeds_in) 

pler_phyt <- phyt %>%
  filter(phytometer == "PLER") %>%
  rename(seeds_out = seeds_mat) %>%
  select(-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes,-biomass_g) %>%
  rename(species = phytometer, background = seed_sp) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_")) %>%
  inner_join(dens)

## LAPL phytometer and and background seed/biomass production
lapl_back <- back %>%
  filter(seed_sp == "LAPL") %>%
  mutate(seeds_out = seeds_mat + seeds_immat, background = "LAPL", PLER_seeds_in = 0, BRHO_seeds_in = 0, FEMI_seeds_in = 0, waterN_treatment = paste(trt_water, trt_N, sep = "_")) %>%
  select(-seeds_mat,-biomass_g,-seeds_immat,-empty_glumes,-seedmass_g,-glumemass_g,-nodes,-senescence) %>%
  rename(species = seed_sp, LAPL_seeds_in = seeds_in)


lapl_phyt <- phyt %>%
  filter(phytometer == "LAPL") %>%
  mutate(seeds_out = seeds_mat + seeds_immat) %>%
  select(-seeds_mat,-biomass_g,-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes) %>%
  rename(species = phytometer, background = seed_sp) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_")) %>%
  inner_join(dens)

## Join all species
brho_dat <- full_join(brho_back,brho_phyt)

brho_dat2 <- full_join(brho_back2,brho_phyt2)

femi_dat <- full_join(femi_back,femi_phyt)

pler_dat <- full_join(pler_back,pler_phyt)

lapl_dat <- full_join(lapl_back,lapl_phyt)

join1 <- full_join(brho_dat,femi_dat)

join1a <- full_join(brho_dat2,femi_dat)

join2 <- full_join(pler_dat,lapl_dat)

model_dat <- full_join(join1,join2)

model_dat2 <- full_join(join1a,join2)



