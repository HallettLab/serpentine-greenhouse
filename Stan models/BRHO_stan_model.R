## Loads rstan and sets the number of useable cores based on your computer
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read in data
data <- read.csv(paste(datpath, "/model_dat2.csv", sep = "")) %>%
  select(-X)

#data <- read.csv("model_dat_0.csv", sep = ",") %>%
#  select(-X)

## Set initial parameter estimates, once for each chain to run
#initials <- list(lambda=10, alpha_pler=0.01, alpha_brho=0.01, alpha_lapl=0.01,
#                 alpha_femi=0.01)


## Subset data for competitor and treatment of interest
dat <- subset(data, species == "BRHO")
dat <- subset(dat, waterN_treatment == "lo_hi")

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))

### Set population context for each species as seeds in
pler <- as.integer(dat$PLER_seeds_in)
brho <- as.integer(dat$BRHO_seeds_in)
lapl <- as.integer(dat$LAPL_seeds_in)
femi <- as.integer(dat$FEMI_seeds_in)

### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- brho

Plot <- dat$block

######################################################

# high high initials
initials <- list(lambda=235.2, alpha_pler=0.22, alpha_brho=1.09, alpha_lapl=0.24,
                 alpha_femi=0.32, epsilon=rep(1,P), sigma = 311.5)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_brho_hi_hi <- stan(file = "Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_brho_hi_hi, pars="lambda")
pairs(no_dist_seeds_brho_hi_hi)

### Save posterior distributions to file
save(no_dist_seeds_brho_hi_hi, file = "no_dist_seeds_brho_hi_hi")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_brho_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
brho_hi_hi <- rstan::extract(no_dist_seeds_brho_hi_hi)
acf(brho_hi_hi$lambda)

######################################################

# high intermediate initials
initials <- list(lambda=52.55, alpha_pler=0.16, alpha_brho=0.83, alpha_lapl=0.19,
                 alpha_femi=0.37, epsilon=rep(1,P), sigma = 172.19)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_brho_hi_int <- stan(file = "Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_brho_hi_int, pars="lambda")
pairs(no_dist_seeds_brho_hi_int)

### Save posterior distributions to file
save(no_dist_seeds_brho_hi_int, file = "no_dist_seeds_brho_hi_int")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_brho_hi_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
brho_hi_int <- rstan::extract(no_dist_seeds_brho_hi_int)
acf(brho_hi_int$lambda)

######################################################

# high low initials
initials <- list(lambda=9.5, alpha_pler=0.1, alpha_brho=0.5, alpha_lapl=0.05,
                 alpha_femi=0.4, epsilon=rep(1,P), sigma = 200)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_brho_hi_lo <- stan(file = "Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_brho_hi_lo, pars="lambda")
pairs(no_dist_seeds_brho_hi_lo)

### Save posterior distributions to file
save(no_dist_seeds_brho_hi_lo, file = "no_dist_seeds_brho_hi_lo")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_brho_hi_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
brho_hi_lo <- rstan::extract(no_dist_seeds_brho_hi_lo)
acf(brho_hi_lo$lambda)

######################################################

# low high initials
initials <- list(lambda=372, alpha_pler=0.5, alpha_brho=3.1, alpha_lapl=0.4,
                 alpha_femi=0.4, epsilon=rep(1,P), sigma = 569)
initials1<- list(initials, initials, initials)

no_dist_seeds_brho_lo_hi <- stan(file = "Rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_brho_lo_hi, pars="lambda")
pairs(no_dist_seeds_brho_lo_hi)

### Save posterior distributions to file
save(no_dist_seeds_brho_lo_hi, file = "no_dist_seeds_brho_lo_hi")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_brho_lo_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
brho_lo_hi <- rstan::extract(no_dist_seeds_brho_lo_hi)
acf(brho_lo_hi$lambda)

######################################################

# low intermediate initials
initials <- list(lambda=25, alpha_pler=0.08, alpha_brho=0.57, alpha_lapl=0.04,
                 alpha_femi=0.07, epsilon=rep(1,P), sigma = 18)
initials1<- list(initials, initials, initials)

no_dist_seeds_brho_lo_int <- stan(file = "Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot"),
                                 iter = 10000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_brho_lo_int, pars="lambda")
pairs(no_dist_seeds_brho_lo_int)

### Save posterior distributions to file
save(no_dist_seeds_brho_lo_int, file = "no_dist_seeds_brho_lo_int")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_brho_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
brho_lo_int <- rstan::extract(no_dist_seeds_brho_lo_int)
acf(brho_lo_int$lambda)

######################################################

# low low initials
initials <- list(lambda=0.28, alpha_pler=0.02, alpha_brho=0.003, alpha_lapl=-0.001,
                 alpha_femi=0.12, epsilon=rep(1,P), sigma = 136)
initials1<- list(initials, initials, initials)

no_dist_seeds_brho_lo_lo <- stan(file = "Constrained_rplot_four_species_BH_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot"),
                                  iter = 12000, chains = 3, thin = 5, control = list(adapt_delta = 0.9999999999999999, max_treedepth =50),
                                  init = initials1)

traceplot(no_dist_seeds_brho_lo_lo, pars="lambda")
pairs(no_dist_seeds_brho_lo_lo)

### Save posterior distributions to file
save(no_dist_seeds_brho_lo_lo, file = "no_dist_seeds_brho_lo_lo")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_brho_lo_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
brho_lo_lo <- rstan::extract(no_dist_seeds_brho_lo_lo)
acf(brho_lo_lo$lambda)