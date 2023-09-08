## Loads rstan and sets the number of useable cores based on your computer
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

bg <- 0.98
pg <- 0.92
lg <- 0.32

## Read in data
data <- read.csv(paste(datpath, "/model_dat2_3spp.csv", sep = "")) %>%
  select(-X) %>%
  mutate(BRHO_seeds_in = ifelse(BRHO_seeds_in %in% 1, 1/0.98, BRHO_seeds_in)) %>%
  mutate(PLER_seeds_in = ifelse(PLER_seeds_in %in% 1, 1/0.92, PLER_seeds_in)) %>%
  mutate(LAPL_seeds_in = ifelse(LAPL_seeds_in %in% 1, 1/0.32, LAPL_seeds_in))

## Subset data for competitor and treatment of interest
dat <- subset(data, species == "BRHO")
dat_alone <- dat %>%
  subset(waterN_treatment == "hi_hi") %>%
  subset(background == "none") %>%
  na.omit()

dat_tog <- dat %>%
  subset(waterN_treatment == "hi_hi") %>%
  subset(background != "none") %>%
  na.omit()

#### Fit model alone ####
## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- dat_alone$seeds_out
N <- as.integer(length(Fecundity_alone))
intra <- dat_alone$BRHO_seeds_in

initials <- list(lambda=100)
initials1 <- list(initials, initials, initials)

brho_hi_hi_alone <- stan(file = "brho_alone_model.stan", 
                   data = c("N", "Fecundity", "intra","bg"),
                   iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =10),
                   init = initials1)

# lambda mean = 112.86, sd = 5.32

f_mean <- 112.86
f_sd <- 5.32

### Set population context for each species as seeds in
pler <- dat_tog$PLER_seeds_in
brho <- dat_tog$BRHO_seeds_in
lapl <- dat_tog$LAPL_seeds_in

Fecundity <- dat_tog$seeds_out
N <- as.integer(length(Fecundity))
intra <- brho

### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))
Plot <- dat_tog$block


######################################################
# HIGH HIGH
# high high dirty run
initials <- list(lambda=113, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(.3,P), sigma = 13)
initials1 <- list(initials, initials, initials)

brho_hi_hi <- stan(file = "brho_community_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg", "f_mean", "f_sd"),
                   iter = 50000, chains = 3, thin = 6, control = list(adapt_delta = 0.99, max_treedepth =50),
                   init = initials1)

### Save posterior distributions to file
save(brho_hi_hi, file = "brho_hi_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(brho_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))
pairs(brho_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
brho_hi_hi <- rstan::extract(brho_hi_hi)
acf(brho_hi_hi$lambda)

cor(brho_hi_hi$alpha_pler, brho_hi_hi$alpha_lapl)
