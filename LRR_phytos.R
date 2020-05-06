###############################
######Log Response Ratio#######
###############################

##load libraries
library(tidyverse)

######PLER seed LRR with competition######
##read in data
seeds_phyt <- read.csv(paste(datpath, "/seed_phytometers.csv", sep = "")) %>%
  filter(!phytometer == "FEMI") %>%
  filter(!potential_seed == "NA")

##pler seed production LRR data manipulation
none_pler <- seeds_phyt %>%
  filter(phytometer == "PLER") %>%
  filter(seed_sp == "none" | seed_density == "none") %>%
  group_by(trt_water, trt_N) %>%
  summarize(nocomp_seed = mean(potential_seed)) %>%
  mutate(nocomp_seed = nocomp_seed+1)

brho_pler <- seeds_phyt %>%
  filter(phytometer == "PLER") %>%
  filter(seed_sp == "BRHO") %>%
  mutate(potential_seed = potential_seed+1)

comp_pler <- seeds_phyt %>%
  filter(phytometer == "PLER") %>%
  filter(!seed_sp == "none") %>%
  mutate(potential_seed = potential_seed+1)

LRR_pler_brho <- left_join(none_pler,brho_pler) %>%
  mutate(LRR = log(potential_seed/nocomp_seed))

LRR_pler <- left_join(none_pler,comp_pler) %>%
  mutate(LRR = log(potential_seed/nocomp_seed))

##plot pler seed production LRR
ggplot(LRR_pler) + geom_boxplot(aes(trt_N,LRR,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_grid(seed_sp~seed_density) +
  ylab("Plantago seed production LRR")

ggplot(LRR_pler_brho) + geom_boxplot(aes(trt_N,LRR,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_wrap(~seed_density) +
  ylab("Plantago seed production LRR")

######BRHO seed LRR with competition######
##brho seed production LRR data manipulation
none_brho <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(seed_sp == "none" | seed_density == "none") %>%
  group_by(trt_water, trt_N) %>%
  summarize(nocomp_matureseed = mean(num_seeds),
            nocomp_allseed = mean(potential_seed)) %>%
  mutate(nocomp_matureseed = nocomp_matureseed+1,
         nocomp_allseed = nocomp_allseed+1)

pler_brho <- seeds_phyt %>%
  filter(phytometer == "BRHO")%>%
  filter(seed_sp == "PLER") %>%
  select(-comment,-biomass_mg) %>%
  mutate(num_seeds = num_seeds+1, potential_seed = potential_seed+1)

comp_brho <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(!seed_sp == "none") %>%
  select(-comment,-biomass_mg) %>%
  mutate(num_seeds = num_seeds+1, potential_seed = potential_seed+1)

LRR_brho_pler <- left_join(none_brho,pler_brho) %>%
  mutate(LRR_mature = log(num_seeds/nocomp_matureseed),
         LRR_all = log(potential_seed/nocomp_allseed))

LRR_brho <- left_join(none_brho,comp_brho) %>%
  mutate(LRR_mature = log(num_seeds/nocomp_matureseed),
         LRR_all = log(potential_seed/nocomp_allseed))

##plot brho seed production LRR
ggplot(LRR_brho_pler) + geom_boxplot(aes(trt_N,LRR_all,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_wrap(~seed_density) +
  ylab("Bromus seed production LRR")

ggplot(LRR_brho_pler) + geom_boxplot(aes(trt_N,LRR_mature,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_wrap(~seed_density) +
  ylab("Bromus seed production LRR")

ggplot(LRR_brho) + geom_boxplot(aes(trt_N,LRR_all,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_grid(seed_sp~seed_density) +
  ylab("Plantago seed production LRR")

ggplot(LRR_brho) + geom_boxplot(aes(trt_N,LRR_mature,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_grid(seed_sp~seed_density) +
  ylab("Plantago seed production LRR")

######LAPL seed LRR with competition######
##lapl seed production LRR data manipulation
none_lapl <- seeds_phyt %>%
  filter(phytometer == "LAPL") %>%
  filter(seed_sp == "none" | seed_density == "none") %>%
  group_by(trt_water, trt_N) %>%
  summarize(nocomp_matureseed = mean(num_seeds),
            nocomp_allseed = mean(potential_seed)) %>%
  mutate(nocomp_matureseed = nocomp_matureseed+1,
         nocomp_allseed = nocomp_allseed+1)

comp_lapl <- seeds_phyt %>%
  filter(phytometer == "LAPL") %>%
  filter(!seed_sp == "none") %>%
  select(-comment,-biomass_mg) %>%
  mutate(num_seeds = num_seeds+1, potential_seed = potential_seed+1)

LRR_lapl <- left_join(none_lapl,comp_lapl) %>%
  mutate(LRR_mature = log(num_seeds/nocomp_matureseed),
         LRR_all = log(potential_seed/nocomp_allseed))

##plot lapl seed production LRR
ggplot(LRR_lapl) + geom_boxplot(aes(trt_N,LRR_all,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_grid(seed_sp~seed_density) +
  ylab("Layia seed production LRR")

ggplot(LRR_lapl) + geom_boxplot(aes(trt_N,LRR_mature,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_grid(seed_sp~seed_density) +
  ylab("Layia seed production LRR")
