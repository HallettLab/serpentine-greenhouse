###############################################
######Fitting biomass and seed production######
###############################################

library(tidyverse)

## Read in data
back <- read.csv(paste(datpath, "/seed_biomass_background.csv", sep = "")) 
phyt <- read.csv(paste(datpath, "/seed_biomass_phytometers.csv", sep = ""))

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
  filter(senescence == "yes") 

ggplot(brho_back_bio) + geom_point(aes(totalmass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(totalmass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

lm_brho_back_bio <- lm(total_seed~totalmass_g,dat=brho_back_bio)
summary(lm_brho_back_bio)
#rsquared = 0.8839
#estimated seeds = (110.33*biomass_g) - 27.14
#taking out 10 samples/pots improved the linear relationship between seed production and biomass

# BRHO samples that did not senesce
brho_back_nosen <- brho_back %>%
  filter(senescence == "no") 
#mostly high N (9/10 samples), but 5 are high water and 4 are lo water, so water is probably not the biggest factor
#5/10 samples removed are from B3, it might be a block issue? 2 from B2 and 3 from B4
#1 int N sample, high water

# Recalculate seed production for brho background based on linear equation
brho_back_adj <- brho_back %>%
  mutate(seeds_out_new = trunc((110.33*totalmass_g)-27.14))
brho_back_adj$seeds_out_new[brho_back_adj$seeds_out_new < 0] <- 0

brho_back_adj <- brho_back_adj %>%
  select(-seeds_mat,-seeds_immat,-senescence,-seedmass_g,-biomass_g,-empty_glumes,-glumemass_g,-nodes) %>%
  rename(seeds_out = total_seed)

write.csv(brho_back_adj, "brho_back2.csv")
  
## BRHO phytometer biomass and seed - all blocks
brho_phyt <- phyt %>%
  filter(phytometer == "BRHO") %>%
  select(-nodes,-comment) %>%
  mutate(total_seed = seeds_mat + seeds_immat)

ggplot(brho_phyt) + geom_point(aes(biomass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(biomass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

brho_phyt_lm <- lm(total_seed~biomass_g,data=brho_phyt)
summary(brho_phyt_lm)
#rsquared really low, very bad linear fit

# After accounting for senescence and presence of reproductive parts
brho_phyt_sen_a <- brho_phyt %>%
  filter(potential_seed>0,empty_glumes>0,total_seed>0)

brho_phyt_sen_b <- brho_phyt %>%
  filter(empty_glumes == 0, potential_seed > 0)

brho_phyt_sen <- full_join(brho_phyt_sen_a,brho_phyt_sen_b) %>%
  filter(pot_id != 97)
#took out outlier pot 97
#the filter above keeps samples that fully developed seed (no empty glumes)
#and samples that developed reproductive parts

ggplot(brho_phyt_sen) + geom_point(aes(biomass_g,total_seed,shape=trt_water,color=trt_N)) +
  geom_smooth(aes(biomass_g,total_seed),method = "lm",se=F) +
  scale_shape_manual(values = c(19,1))

lm_brho_phyt_sen <- lm(total_seed~biomass_g,dat=brho_phyt_sen)
summary(lm_brho_phyt_sen)
#rsquared = 0.8898
#estimated seeds = (180.7284*biomass_g) - 0.2494

# Recalculate seed production for brho phytometers based on linear equation
brho_phyt_adj <- brho_phyt %>%
  mutate(seeds_out_new = trunc((180.7284*biomass_g)-0.2494), seeds_in=1) %>%
  select(-seeds_mat,-seeds_immat,-potential_seed,-empty_glumes,) %>%
  rename(seeds_out = total_seed)

write.csv(brho_phyt_adj, "brho_phyt2.csv")

