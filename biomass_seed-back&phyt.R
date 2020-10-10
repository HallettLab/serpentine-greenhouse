###############################################
######Fitting biomass and seed production######
###############################################

library(tidyverse)

## Read in data
back <- read.csv(paste(datpath, "/seed_biomass_background.csv", sep = "")) 
phyt <- read.csv(paste(datpath, "/seed_biomass_phytometers.csv", sep = ""))

## FEMI background biomass and seed - Block1
femi_back_1 <- back %>%
  filter(seed_sp == "FEMI") %>%
  filter(block == 1) %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  mutate(out_in_nodes = nodes/seeds_in) %>%
  mutate(out_in_seeds = total_seed/seeds_in) %>%
  mutate(total_mass_g = biomass_g + seedmass_g + glumemass_g) %>%
  select(-biomass_g,-seedmass_g,-glumemass_g,-empty_glumes) 
# Check out seeds vs nodes
ggplot(femi_back_1) + geom_point(aes(trt_N,out_in_nodes,fill=trt_water)) +
  facet_wrap(~seed_density)
ggplot(femi_back_1) + geom_point(aes(trt_N,out_in_seeds,fill=trt_water)) +
  facet_wrap(~seed_density)

ggplot(femi_back_1) + geom_point(aes(total_mass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(total_mass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))
ggplot(femi_back_1) + geom_point(aes(total_mass_g,nodes,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(total_mass_g,nodes),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

lm_femiback1_seed <- lm(total_seed~total_mass_g, dat=femi_back_1)
#rsquared = 0.91
lm_femiback1_nodes <- lm(nodes~total_mass_g,dat=femi_back_1)
#rsquared = 0.962
#intercept = 16.5; coefficient total_mass_g = 189.9
# Conclusion: use nodes with the following linear equation:
# nodes(seed) = (189.9*biomass_g) + 16.5
# Note: seed fell during the experiment and sample collection process, 
#       therefore, # of nodes are more reliable and also a better fit

## FEMI background seed production calculated with biomass - Blocks 2,3,4
femi_back_234 <- back %>%
  filter(seed_sp == "FEMI") %>%
  filter(block != 1) %>%
  mutate(nodes = ((biomass_g*189.9)+16.5)) %>%
  select(-seedmass_g,-glumemass_g,-seeds_mat,-seeds_immat,-empty_glumes) %>%
  


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
# After accounting for senescence
brho_back_sen <- read.csv(paste(datpath, "/bromus_back_phen.csv", sep = ""))
# Biomass senescence
brho_back_bio <- brho_back_sen %>%
  filter(senes_biomass == "yes") %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  mutate(totalmass_g = biomass_g + seedmass_g)

ggplot(brho_back_bio) + geom_point(aes(totalmass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(totalmass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

lm_brho_back_bio <- lm(total_seed~totalmass_g,dat=brho_back_bio)
# estimated seeds = (110.33*biomass_g) - 27.14
# taking out 10 samples/pots improved the linear relationship between seed production and biomass
# data frame containing samples that were removed for best fit line
brho_back_nosen <- brho_back_sen %>%
  filter(senes_biomass == "no") 
# mostly high N (9/10 samples), but 5 are high water and 4 are lo water, so water is probably not the biggest factor
# 5/10 samples removed are from B3, it might be a block issue? 2 from B2 and 3 from B4
# 1 int N sample, high water
# Recalculate seed production based on linear equation
brho_back_nosen <- brho_back_sen %>%
  filter(senes_biomass == "no") %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  mutate(totalmass_g = biomass_g + seedmass_g) %>%
  mutate(total_seed_new = (110.33*totalmass_g)-27.14)


## BRHO phytometer biomass and seed - all blocks
brho_phyt <- phyt %>%
  filter(phytometer == "BRHO")
brho_phyt$immat_seed[is.na(brho_phyt$immat_seed)] = 0

brho_phyt <- brho_phyt %>%
  select(-glumes,-comment) %>%
  mutate(total_seed = mat_seed + immat_seed)

ggplot(brho_phyt) + geom_point(aes(biomass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(biomass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

brho_phyt_lm <- lm(total_seed~biomass_g,data=brho_phyt)
summary(brho_phyt_lm)
# After accounting for senescence
brho_phyt_sen <- read.csv(paste(datpath, "/bromus_phytos_phen.csv", sep = "")) %>%
  filter(senescence == "yes") %>%
  mutate(total_seed = mat_seed + immat_seed)

ggplot(brho_phyt_sen) + geom_point(aes(biomass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(biomass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

lm_brho_back_bio <- lm(total_seed~totalmass_g,dat=brho_back_bio)
# 27% (46/168 phytometers) of bromus phytometers senesced 50% before end of experiment
# Accounting for senescence does not improve seed output as a function of biomass
# But I counted all the seed, so no need to incorporate biomass

