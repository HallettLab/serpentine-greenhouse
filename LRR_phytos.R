###############################
######Log Response Ratio#######
###############################

##load libraries
library(tidyverse)

######PLER seed LRR with BRHO competition######
##read in data
seeds_phyt <- read.csv(paste(datpath, "/seed_phytometers.csv", sep = "")) %>%
  filter(!phytometer == "FEMI") %>%
  filter(!potential_seed == "NA")

##pler seed production LRR data manipulation
none_pler <- seeds_phyt %>%
  filter(phytometer == "PLER") %>%
  filter(seed_sp == "none" | seed_density == "none") %>%
  group_by(trt_water, trt_N) %>%
  summarize(nocomp_seed = mean(num_seeds)) 

brho_pler <- seeds_phyt %>%
  filter(phytometer == "PLER")%>%
  filter(seed_sp == "BRHO") %>%
  select(-potential_seed,-comment,-biomass_mg)

LRR_pler <- left_join(none_pler,brho_pler) %>%
  mutate(LRR = log(num_seeds/nocomp_seed))

##plot pler seed production LRR
ggplot(LRR_pler) + geom_boxplot(aes(trt_N,LRR,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_wrap(~seed_density) +
  ylab("Plantago seed production LRR")


######BRHO seed LRR with PLER competition######
##brho seed production LRR data manipulation
none_brho <- seeds_phyt %>%
  filter(phytometer == "BRHO") %>%
  filter(seed_sp == "none" | seed_density == "none") %>%
  group_by(trt_water, trt_N) %>%
  summarize(nocomp_matureseed = mean(num_seeds),
            nocomp_allseed = mean(potential_seed)) 

pler_brho <- seeds_phyt %>%
  filter(phytometer == "BRHO")%>%
  filter(seed_sp == "PLER") %>%
  select(-comment,-biomass_mg)

LRR_brho <- left_join(none_brho,pler_brho) %>%
  mutate(LRR_mature = log(num_seeds/nocomp_matureseed),
         LRR_all = log(potential_seed/nocomp_allseed))

##plot brho seed production LRR
ggplot(LRR_brho) + geom_boxplot(aes(trt_N,LRR_all,fill=trt_water)) + 
  geom_hline(yintercept=0,linetype="dashed") + facet_wrap(~seed_density) +
  ylab("Bromus seed production LRR")
