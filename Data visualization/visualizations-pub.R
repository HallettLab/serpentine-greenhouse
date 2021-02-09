####################################
##### Publication ready graphs #####
####################################

library(tidyverse)

## Read in data
seed_biomass_dat <- read.csv(paste(datpath, "/clean_dat2.csv", sep = ""))
stems_dat <-  read.csv(paste(datpath, "/stems_background.csv", sep = ""))

## SE function
calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

## Set theme
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 12),
              strip.text= element_text(size = 12),
              axis.text = element_text(size = 12))

########################
## Background biomass ##
########################

## Data manipulation 
pler_back <- seed_biomass_dat %>% 
  filter(species == "PLER", background == "PLER")

femi_back <- seed_biomass_dat %>%
  filter(species == "FEMI", background == "FEMI")

lapl_back <- seed_biomass_dat %>%
  filter(species == "LAPL", background == "LAPL")

brho_back <- seed_biomass_dat %>%
  filter(species == "BRHO", background == "BRHO")

femi_pler_back <- full_join(pler_back,femi_back)
brho_lapl_back <- full_join(lapl_back,brho_back)

back_biomass <- full_join(femi_pler_back,brho_lapl_back) %>%
  select(-X,-seeds_in,-seeds_out,-out_in)
  
## Reorder and rename factors
back_biomass$trt_water <- factor(back_biomass$trt_water , levels = c("lo","hi"))
back_biomass$trt_N <- factor(back_biomass$trt_N , levels = c("lo","int","hi"))

water.labs <- c("Dry", "Wet")
names(water.labs) <- c("lo", "hi")
n.labs <- c("Low","Intermediate","High")
names(n.labs) <- c("lo", "hi")
sp.labs <- c("Bromus hordeaceus", "Festuca microstachys","Layia platyglossa","Plantago erecta")
names(sp.labs) <- c("BRHO", "FEMI","LAPL","PLER")
dens.labs <- c("Low seed addition", "High seed addition")
names(dens.labs) <- c("lo", "hi")

## Data visualization
ggplot(back_biomass) + geom_boxplot(aes(trt_N,biomass_g,fill=trt_water)) + 
  facet_grid(seed_density~background,
             labeller = labeller(background=sp.labs,seed_density=dens.labs)) +
  ylab("Biomass (g)") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Water treatments", labels = c("Dry","Wet"), 
                    values=c("grey80","grey40")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(strip.text.x = element_text(face = "italic"))
  
#################################
## Background stem/recruitment ##
#################################

## Data manipulation
rec_dat <- stems_dat %>%
  filter(seed_sp != "none") %>%
  mutate(recruitment = (stem_density/seed_added)*100) %>%
  rename(background = seed_sp)

## Reorder and rename factors
rec_dat$trt_water <- factor(rec_dat$trt_water , levels = c("lo","hi"))
rec_dat$trt_N <- factor(rec_dat$trt_N , levels = c("lo","int","hi"))

## Data visualization
# Recruitment
ggplot(rec_dat) + geom_boxplot(aes(trt_N,recruitment,fill=trt_water)) + 
  facet_grid(seed_density~background,
             labeller = labeller(background=sp.labs,seed_density=dens.labs)) +
  ylab("Recruitment (%)") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Water treatments", labels = c("Dry","Wet"), 
                    values=c("grey80","grey40")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(strip.text.x = element_text(face = "italic"))

# Stem density
ggplot(rec_dat) + geom_boxplot(aes(trt_N,stem_density,fill=trt_water)) + 
  facet_grid(seed_density~background,
             labeller = labeller(background=sp.labs,seed_density=dens.labs), 
             scale = "free") +
  ylab("Stem density") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Water treatments", labels = c("Dry","Wet"), 
                    values=c("grey80","grey40")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(strip.text.x = element_text(face = "italic"))

