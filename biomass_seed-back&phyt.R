###############################################
######Fitting biomass and seed production######
###############################################

library(tidyverse)

## Read in data
back <- read.csv(paste(datpath, "/seed_biomass_background.csv", sep = "")) 
phyt <- read.csv(paste(datpath, "/seed_biomass_phytometers.csv", sep = ""))

## FEMI seed production based on biomass ##
# FEMI background biomass and seed - Block1
femi_back_1 <- back %>%
  filter(seed_sp == "FEMI") %>%
  filter(block == 1) %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  mutate(out_in_nodes = nodes/seeds_in) %>%
  mutate(out_in_seeds = total_seed/seeds_in) %>%
  mutate(total_mass_g = biomass_g + seedmass_g + glumemass_g) %>%
  select(-biomass_g,-seedmass_g,-glumemass_g,-empty_glumes) 
#Check out seeds vs nodes
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
#Conclusion: estimate nodes with the following linear equation:
#nodes(seed) = (189.9*biomass_g) + 16.5
#Note: seed fell during the experiment and sample collection process, 
#       therefore, # of nodes are more reliable and also a better fit

# FEMI background seed production calculated with biomass - Blocks 2,3,4
femi_back_234 <- back %>%
  filter(seed_sp == "FEMI") %>%
  filter(block != 1) %>%
  mutate(nodes = ((biomass_g*189.9)+16.5)) %>%
  select(-seedmass_g,-glumemass_g,-seeds_mat,-seeds_immat,-empty_glumes) %>%
  
# Join Block1 and Blocks 2,3,4 of FEMI background seed production
femi_back1 <- femi_back_1 %>%
  rename(biomass_g = total_mass_g)
#clean up
femi_back_clean <- full_join(femi_back1,femi_back_234) %>%
  select(-seeds_mat,-seeds_immat,-total_seed,-out_in_seeds) %>%
  rename(seeds_out = nodes, out_in = out_in_nodes) %>%
  mutate(out_in = seeds_out/seeds_in)


## BRHO background seed production based on senescence ##
# BRHO seed production vs. biomass
#Look at data before filtering for senescence
brho_back <- back %>%
  filter(seed_sp == "BRHO") %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  mutate(totalmass_g = biomass_g + seedmass_g)

ggplot(brho_back) + geom_boxplot(aes(trt_N,totalmass_g,fill=trt_water)) +
  facet_wrap(~seed_density)

ggplot(brho_back) + geom_point(aes(totalmass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(totalmass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

lm_brho_back <- lm(total_seed~totalmass_g,dat=brho_back)
summary(lm_brho_back)
#rsquared = 0.426

# BRHO samples that senesced
#Note: went back to samples to look at greenness
#       samples that had more than 50% senescence are "yes" for senescence
brho_back_bio <- brho_back %>%
  filter(senescence == "yes") %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  mutate(totalmass_g = biomass_g + seedmass_g)

ggplot(brho_back_bio) + geom_point(aes(totalmass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(totalmass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

lm_brho_back_bio <- lm(total_seed~totalmass_g,dat=brho_back_bio)
summary(lm_brho_back_bio)
#rsquared = .8839
#estimated seeds = (110.33*biomass_g) - 27.14
#taking out 10 samples/pots improved the linear relationship between seed production and biomass

# BRHO samples that did not senesce
brho_back_nosen <- brho_back %>%
  filter(senescence == "no") 
#mostly high N (9/10 samples), but 5 are high water and 4 are lo water, so water is probably not the biggest factor
#5/10 samples removed are from B3, it might be a block issue? 2 from B2 and 3 from B4
#1 int N sample, high water

# Recalculate seed production for based on linear equation
brho_back_nosen <- brho_back_nosen %>%
  mutate(total_seed = seeds_mat + seeds_immat) %>%
  mutate(totalmass_g = biomass_g + seedmass_g) %>%
  mutate(total_seed_new = (110.33*totalmass_g)-27.14)

# Join not senesced samples with senesced
#Cleanup
brho_back_sen <- brho_back %>%
  filter(!pot_id %in% c(86,93,140,142,145,149,161,192,207,210)) %>%
  select(-seeds_mat,-seeds_immat,-senescence,-seedmass_g,-biomass_g,-empty_glumes,-glumemass_g,-nodes) %>%
  rename(biomass_g = totalmass_g, seeds_out = total_seed)

brho_back_nosen2 <- brho_back_nosen %>%
  select(-seeds_mat,-seeds_immat,-senescence,-seedmass_g,-biomass_g,-empty_glumes,-total_seed,-glumemass_g,-nodes) %>%
  rename(biomass_g = totalmass_g, seeds_out = total_seed_new)
#Join
brho_back_clean <- full_join(brho_back_sen,brho_back_nosen2)

  
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

