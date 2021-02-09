########################################################
##### Data clean up for visualization and analysis #####
########################################################

library(tidyverse)

source("fitting_seed_biomass.R")

## Read in other data
# seeds_in for phytometers is 0
dens1 <- read.csv(paste(datpath, "/seed_density.csv", sep = ""))
# seeds_in for phytometers is 1
dens2 <- read.csv(paste(datpath, "/seed_density2.csv", sep = ""))

################################
## Clean up for visualization ##
################################

# BRHO phytometer and background seed/biomass production
brho_back_adj1 <- brho_back_adj1 %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"))

brho_phyt_adj1 <- brho_phyt_adj1 %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"))

brho_back_adj2 <- brho_back_adj2 %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"))

brho_phyt_adj2 <- brho_phyt_adj2 %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_")) 

# FEMI phytometer and background seed/biomass production
femi_back_adj <- femi_back_adj %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"), background = "FEMI")

femi_phyt <- phyt %>%
  filter(phytometer == "FEMI") %>%
  rename(seeds_out = seeds_mat) %>%
  select(-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes) %>%
  rename(species = phytometer, background = seed_sp) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"), seeds_in = 1) 

# PLER phytometer and background seed/biomass production
pler_back <- back %>%
  filter(seed_sp == "PLER") %>%
  rename(seeds_out = seeds_mat) %>%
  select(-seeds_immat,-empty_glumes,-seedmass_g,-glumemass_g,-nodes,-senescence) %>%
  rename(species = seed_sp) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"), background = "PLER")

pler_phyt <- phyt %>%
  filter(phytometer == "PLER") %>%
  rename(seeds_out = seeds_mat) %>%
  select(-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes) %>%
  rename(species = phytometer, background = seed_sp) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"), seeds_in = 1) 

# LAPL phytometer and background seed/biomass production
lapl_back <- back %>%
  filter(seed_sp == "LAPL") %>%
  mutate(seeds_out = seeds_mat + seeds_immat, background = "LAPL") %>%
  select(-seeds_mat,-seeds_immat,-empty_glumes,-seedmass_g,-glumemass_g,-nodes,-senescence) %>%
  rename(species = seed_sp) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_")) 

lapl_phyt <- phyt %>%
  filter(phytometer == "LAPL") %>%
  mutate(seeds_out = seeds_mat + seeds_immat) %>%
  select(-seeds_mat,-seeds_immat,-potential_seed,-empty_glumes,-comment,-nodes) %>%
  rename(species = phytometer, background = seed_sp) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"), seeds_in=1)

# Join all species
brho_dat1 <- full_join(brho_back_adj1,brho_phyt_adj1)
brho_dat2 <- full_join(brho_back_adj2,brho_phyt_adj2)
femi_dat <- full_join(femi_back_adj,femi_phyt)
pler_dat <- full_join(pler_back,pler_phyt)
lapl_dat <- full_join(lapl_back,lapl_phyt)
join1 <- full_join(brho_dat1,femi_dat)
join1a <- full_join(brho_dat2,femi_dat)
join2 <- full_join(pler_dat,lapl_dat)

# Clean_dat1 consists of bromus phyt and back "not senensced" corrected samples
clean_dat1 <- full_join(join1,join2) %>%
  mutate(out_in = seeds_out/seeds_in)

# Clean_dat2 consists of bromus phyt and back all corrected samples
clean_dat2 <- full_join(join1a,join2) %>%
  mutate(out_in = seeds_out/seeds_in)

# CSV file
#write.csv (clean_dat1, "clean_dat1.csv")
#write.csv(clean_dat1, "clean_dat2.csv")

## Clean up for models with seeds in as zero
#BRHO phytometer and background seed/biomass production
brho_back_adj1 <- brho_back_adj1 %>%
  mutate(PLER_seeds_in = 0, LAPL_seeds_in = 0, FEMI_seeds_in = 0) %>%
  rename(BRHO_seeds_in = seeds_in)

brho_phyt_adj1 <- brho_phyt_adj1 %>%
  inner_join(dens1) %>%
  select(-seeds_in)

brho_back_adj2 <- brho_back_adj2 %>%
  mutate(PLER_seeds_in = 0, LAPL_seeds_in = 0, FEMI_seeds_in = 0) %>%
  rename(BRHO_seeds_in = seeds_in)

brho_phyt_adj2 <- brho_phyt_adj2 %>%
  inner_join(dens1) %>%
  select(-seeds_in)

# FEMI phytometer and background seed/biomass production
femi_back_adj <- femi_back_adj %>%
  mutate(PLER_seeds_in = 0, LAPL_seeds_in = 0, BRHO_seeds_in = 0) %>%
  rename(FEMI_seeds_in = seeds_in)

femi_phyt <- femi_phyt %>%
  inner_join(dens1) %>%
  select(-seeds_in)

# PLER phytometer and background seed/biomass production
pler_back <- pler_back %>%
  mutate(BRHO_seeds_in = 0, LAPL_seeds_in = 0, FEMI_seeds_in = 0) %>%
  rename(PLER_seeds_in = seeds_in)

pler_phyt <- pler_phyt %>%
  inner_join(dens1) %>%
  select(-seeds_in)

# LAPL phytometer and and background seed/biomass production
lapl_back <- lapl_back %>%
  mutate(PLER_seeds_in = 0, BRHO_seeds_in = 0, FEMI_seeds_in = 0) %>%
  rename(LAPL_seeds_in = seeds_in)

lapl_phyt <- lapl_phyt %>%
  inner_join(dens1) %>%
  select(-seeds_in)

# Join all species
brho_dat <- full_join(brho_back_adj1,brho_phyt_adj1)
brho_dat2 <- full_join(brho_back_adj2,brho_phyt_adj2)
femi_dat <- full_join(femi_back_adj,femi_phyt)
pler_dat <- full_join(pler_back,pler_phyt)
lapl_dat <- full_join(lapl_back,lapl_phyt)
join1 <- full_join(brho_dat,femi_dat)
join1a <- full_join(brho_dat2,femi_dat)
join2 <- full_join(pler_dat,lapl_dat)
model_dat1 <- full_join(join1,join2) %>%
  select(-biomass_g)
model_dat2 <- full_join(join1a,join2) %>%
  select(-biomass_g)

# Write file
# Some bromus samples adjusted
 #write.csv(model_dat1, "model_dat_1.csv")
# All bromus samples adjusted
 #write.csv(model_dat2, "model_dat_0.csv")

## Clean up for models with seeds in as 1
#BRHO phytometer and background seed/biomass production
#brho_back_adj1 <- brho_back_adj1 %>%
#  mutate(PLER_seeds_in = 1, LAPL_seeds_in = 1, FEMI_seeds_in = 1) %>%
#  rename(BRHO_seeds_in = seeds_in)

#brho_phyt_adj1 <- brho_phyt_adj1 %>%
#  inner_join(dens2) %>%
#  select(-seeds_in)

#brho_back_adj2 <- brho_back_adj2 %>%
#  mutate(PLER_seeds_in = 1, LAPL_seeds_in = 1, FEMI_seeds_in = 1) %>%
#  rename(BRHO_seeds_in = seeds_in)

#brho_phyt_adj2 <- brho_phyt_adj2 %>%
#  inner_join(dens2) %>%
#  select(-seeds_in)

# FEMI phytometer and background seed/biomass production
#femi_back_adj <- femi_back_adj %>%
#  mutate(PLER_seeds_in = 1, LAPL_seeds_in = 1, BRHO_seeds_in = 1) %>%
#  rename(FEMI_seeds_in = seeds_in)

#femi_phyt <- femi_phyt %>%
#  inner_join(dens2) %>%
#  select(-seeds_in)

# PLER phytometer and background seed/biomass production
#pler_back <- pler_back %>%
#  mutate(BRHO_seeds_in = 1, LAPL_seeds_in = 1, FEMI_seeds_in = 1) %>%
#  rename(PLER_seeds_in = seeds_in)

#pler_phyt <- pler_phyt %>%
#  inner_join(dens2) %>%
#  select(-seeds_in)

# LAPL phytometer and and background seed/biomass production
#lapl_back <- lapl_back %>%
#  mutate(PLER_seeds_in = 1, BRHO_seeds_in = 1, FEMI_seeds_in = 1) %>%
#  rename(LAPL_seeds_in = seeds_in)

#lapl_phyt <- lapl_phyt %>%
#  inner_join(dens2) %>%
#  select(-seeds_in)

# Join all species
#brho_dat <- full_join(brho_back_adj1,brho_phyt_adj1)
#brho_dat2 <- full_join(brho_back_adj2,brho_phyt_adj2)
#femi_dat <- full_join(femi_back_adj,femi_phyt)
#pler_dat <- full_join(pler_back,pler_phyt)
#lapl_dat <- full_join(lapl_back,lapl_phyt)
#join1 <- full_join(brho_dat,femi_dat)
#join1a <- full_join(brho_dat2,femi_dat)
#join2 <- full_join(pler_dat,lapl_dat)
#model_dat1 <- full_join(join1,join2) %>%
#  select(-biomass_g)
#model_dat2 <- full_join(join1a,join2) %>%
#  select(-biomass_g)

# Write file
# Some bromus samples adjusted
#write.csv(model_dat1, "model_dat_1.csv")
# All bromus samples adjusted
#write.csv(model_dat2, "model_dat_0.csv")







