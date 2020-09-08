###############################################
######Fitting biomass and seed production######
###############################################

library(tidyverse)

## Read in data
back <- read.csv(paste(datpath, "/seed_biomass_background.csv", sep = "")) 
phyt <- read.csv(paste(datpath, "/seed_biomass_phytometers.csv", sep = ""))


## FEMI background biomass and seed for block 1
femi_back <- back %>%
  filter(seed_sp == "FEMI") %>%
  select(-seedmass_g) %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  filter(block == 1) %>%
  rename(nodes = empty_glumes) %>%
  mutate(out_in = nodes/seeds_in)

ggplot(femi_back) + geom_point(aes(biomass_g,nodes,shape=trt_water,color=trt_N)) +
  scale_shape_manual(values = c(19,1)) + geom_smooth(aes(biomass_g,nodes),method="lm",se=F)

ggplot(femi_back) + geom_point(aes(trt_N,biomass_g,color = trt_water)) +facet_wrap(~seed_density)

ggplot(femi_back) + geom_point(aes(trt_N,out_in,color = trt_water)) +facet_wrap(~seed_density)

femi_lm_nodes <- lm(nodes~biomass_g,data=femi_back)
summary(femi_lm_nodes)
# estimated seeds (nodes) = (16.6*biomass_g) + 236.2


## FEMI phytometer biomass and seed - 2 blocks
femi_phyt <- phyt %>%
  filter(phytometer == "FEMI") %>%
  rename(seed = potential_seed)

ggplot(femi_phyt) + geom_point(aes(biomass_g,seed,color=trt_N,shape=trt_water))+
  scale_shape_manual(values = c(19,1)) + geom_smooth(aes(biomass_g,seed),method="lm",se=F)

femi_phyt_lm <- lm(seed~biomass_g,data=femi_phyt)
summary(femi_phyt_lm)


## BRHO background biomass and seed - all blocks
brho_back <- back %>%
  filter(seed_sp == "BRHO") %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  mutate(totalmass_g = biomass_g + seedmass_g)
  
ggplot(brho_back) + geom_boxplot(aes(trt_N,totalmass_g,fill=trt_water)) +
  facet_wrap(~seed_density)

ggplot(brho_back) + geom_point(aes(totalmass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(totalmass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

## BRHO phytometer biomass and seed - all blocks
brho_phyt <- phyt %>%
  filter(phytometer == "BRHO")
brho_phyt$immat_seed[is.na(brho_phyt$immat_seed)] = 0

brho_phyt <- brho_phyt %>%
  select(-glumes,-comment) %>%
  mutate(total_seed = mat_seed + immat_seed)

ggplot(brho_phyt) + geom_point(aes(biomass_g,total_seed,shape=trt_water,color=trt_N)) +
  scale_shape_manual(values = c(19,1))

brho_phyt_lm <- lm(total_seed~biomass_g,data=brho_phyt)
summary(brho_phyt_lm)



