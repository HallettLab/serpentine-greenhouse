library(tidyverse)
library(minpack.lm)
library(nlstools)

##read in data
cleandat <- read_csv(paste(datpath, "/clean_seed_dat.csv", sep = ""))

seeddat <- cleandat %>%
  select(block:seeds_out, background_comp) %>%
  pivot_longer(seeds_in:seeds_out, "measure", "value") %>%
  mutate(measure2 = paste(individual, measure, sep = "_")) %>%
  select(-measure, -individual) %>%
  pivot_wider(names_from = measure2, values_from = value) %>%
  mutate(waterN_treatment = paste(trt_water, trt_N, sep = "_"))

plerbrho <- seeddat %>%
  filter(background_comp%in%c("BRHO", "PLER", "none"))

# hack to play around with other combos
# plerbrho <- seeddat %>%
#   filter(background_comp%in%c("FEMI", "BRHO", "none")) %>%
#   select(block:background_comp, BRHO_seeds_in, BRHO_seeds_out,FEMI_seeds_in, FEMI_seeds_out, waterN_treatment) %>% 
#   rename(PLER_seeds_in =FEMI_seeds_in,
#          PLER_seeds_out = FEMI_seeds_out)
  


## Set germination and survival fractions from the literature
ps <- .6 # lauren guess
pg <- .92 # gulmon
bs <- .2 # lauren guess
bg <- .98 # gulmon

##### PLANTAGO MODEL ####


# Here, the log of Plantago seeds produced (PLER_seeds_out) is a product of:
# 1) The number of seeds added in relation to 
# 2) The intrinsic growth rate (lambda)
# 3) Competition from fellow Plantago individuals (competition term: aiP; the number of Plantago seeds added is multiplied by the germination fraction, pg)
# 4) Competition from Bromus individuals (competition term: aiB; the number of Bromus seeds added is multiplied by the germination fraction, bg)
m1P <- as.formula(log(PLER_seeds_out +1) ~  log(pg*(PLER_seeds_in+1)*exp(log(lambda)-log((1+aiP*(PLER_seeds_in+1)*pg+aiB*(BRHO_seeds_in+1)*bg)))))


treatments <- unique(plerbrho$waterN_treatment)

Poutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(Poutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1P, start=list(lambda=1, aiP = .01, aiB=.01),
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(plerbrho, waterN_treatment == treatments[i]))
  
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Plantago"
  Poutput <- rbind(Poutput, outreport)
}


#### BROMUS MODEL ###
## With seed bank 

# new model, log transformed
m1B <- as.formula(log(BRHO_seeds_out +1) ~  log(bg*(BRHO_seeds_in+1)*exp(log(lambda)-log((1+aiP*(PLER_seeds_in+1)*pg+aiB*(BRHO_seeds_in+1)*bg)))))


# old model (same model but not log transformed)
#m1 <- as.formula(log(AVseedout +1) ~  log(AVseedin*ag +1)*((lambda)/(1+aiE*log(ERseedin*eg + 1) + aiA*log(AVseedin*ag + 1))))

treatments <- unique(togdat$treatment)

Boutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(Boutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  
  m1out <- nlsLM(m1B, start=list(lambda=1, aiP = .01, aiB=.01), 
                 control=nls.lm.control(maxiter=500), 
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 trace=T,
                 data = subset(plerbrho, waterN_treatment == treatments[i]))
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Bromus"
  Boutput <- rbind(Boutput, outreport)
}

## PUT THE TWO OUTPUTS TOGETHER
model.dat <- rbind(Poutput, Boutput) %>%
  tbl_df() %>%
  dplyr::select(estimate, params, treatment, species) %>%
  spread(params, estimate)

model.dat
 # ggplot(model.dat, aes(x=treatment, y=aiB)) + geom_bar(stat="identity") + facet_wrap(~species)
 # ggplot(model.dat, aes(x=treatment, y=aiP)) + geom_bar(stat="identity") + facet_wrap(~species)


parameter_table <- rbind(Poutput, Boutput) %>%
  tbl_df() %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(estimateout = paste(estimate, "±", se)) %>%
  dplyr::select(estimateout, params, treatment, species) %>%
  spread(params, estimateout) %>%
  select(treatment, species, lambda, aiP, aiB) %>%
  arrange(species, treatment)
