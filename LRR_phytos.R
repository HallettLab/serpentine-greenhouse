###############################
######Log Response Ratio#######
###############################

##load libraries
library(tidyverse)

##read in data
seeds_phyt <- read.csv(paste(datpath, "/seed_phytometers.csv", sep = "")) %>%
  filter(!phytometer == "FEMI") %>%
  filter(!potential_seed == "NA") %>%
  select(-comment,-biomass_mg)
seeds_phyt[seeds_phyt == 0] <- 0.01

##no competition phytometer seed production LRR data manipulation
no_comp <- seeds_phyt %>%
  filter(seed_sp == "none") %>%
  group_by(block,trt_water,trt_N,phytometer) %>%
  summarize(nocomp_seed = potential_seed) 

##competition phytometer seed production LRR data manipulation
comp <- seeds_phyt %>%
  filter(!seed_sp == "none") %>%
  group_by(block,trt_water, trt_N,seed_sp,seed_density,phytometer) %>%
  summarize(comp_seed = potential_seed)

#join datasets for final LRR
LRR <- full_join(no_comp,comp) %>%
  mutate(LRR = log(comp_seed/nocomp_seed))

##Plantago phytometer seed production LRR
LRR_pler <- LRR %>%
  filter(phytometer=="PLER")

ggplot(LRR_pler) + geom_boxplot(aes(trt_N,LRR,fill=trt_water)) +
  facet_grid(seed_sp~seed_density,scale="free") + 
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Plantago phytometer seed production LRR")

##Layia phytometer seed production LRR
LRR_lapl <- LRR %>%
  filter(phytometer=="LAPL")

ggplot(LRR_lapl) + geom_boxplot(aes(trt_N,LRR,fill=trt_water)) +
  facet_grid(seed_sp~seed_density,scale="free") + 
  geom_hline(yintercept=0,linetype="dashed")+
  ylab("Layia phytometer seed production LRR")

##Bromus phytometer seed production LRR
LRR_brho <- LRR %>%
  filter(phytometer=="BRHO")

ggplot(LRR_brho) + geom_boxplot(aes(trt_N,LRR,fill=trt_water)) +
  facet_grid(seed_sp~seed_density,scale="free") + 
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Bromus phytometer seed production LRR")
