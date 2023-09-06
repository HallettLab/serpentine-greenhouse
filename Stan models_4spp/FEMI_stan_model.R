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
dat <- subset(data, species == "FEMI")
dat <- subset(dat, waterN_treatment == "lo_lo")

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
intra <- femi

Plot <- dat$block

fg <- 0.83

######################################################

# high high initials
initials <- list(lambda=178, alpha_pler=0.1, alpha_brho=0.35, alpha_lapl=0.09,
                 alpha_femi=0.31, epsilon=rep(1,P), sigma = 173)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_femi_hi_hi <- stan(file = "FEMI_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","fg"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_femi_hi_hi, pars="lambda")
pairs(no_dist_seeds_femi_hi_hi)

### Save posterior distributions to file
save(no_dist_seeds_femi_hi_hi, file = "no_dist_seeds_femi_hi_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_femi_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
femi_hi_hi <- rstan::extract(no_dist_seeds_femi_hi_hi)
acf(femi_hi_hi$lambda)

######################################################

# high intermediate initials
initials <- list(lambda=13, alpha_pler=0.01, alpha_brho=0.07, alpha_lapl=0.01,
                 alpha_femi=0.06, epsilon=rep(1,P), sigma =119)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_femi_hi_int <- stan(file = "FEMI_Constrained_rplot_four_species_BH_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","fg"),
                                  iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99999999999999, max_treedepth =50),
                                  init = initials1)

traceplot(no_dist_seeds_femi_hi_int, pars="lambda")
pairs(no_dist_seeds_femi_hi_int)

### Save posterior distributions to file
save(no_dist_seeds_femi_hi_int, file = "no_dist_seeds_femi_hi_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_femi_hi_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
femi_hi_int <- rstan::extract(no_dist_seeds_femi_hi_int)
acf(femi_hi_int$lambda)

######################################################

# high low initials
initials <- list(lambda=2, alpha_pler=0.001, alpha_brho=0.02, alpha_lapl=0.002,
                 alpha_femi=0.009, epsilon=rep(1,P), sigma = 181)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_femi_hi_lo <- stan(file = "FEMI_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","fg"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_femi_hi_lo, pars="lambda")
pairs(no_dist_seeds_femi_hi_lo)

### Save posterior distributions to file
save(no_dist_seeds_femi_hi_lo, file = "no_dist_seeds_femi_hi_lo.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_femi_hi_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
femi_hi_lo <- rstan::extract(no_dist_seeds_femi_hi_lo)
acf(femi_hi_lo$lambda)

######################################################

# low high initials
initials <- list(lambda=64, alpha_pler=0.03, alpha_brho=0.16, alpha_lapl=0.02,
                 alpha_femi=0.24, epsilon=rep(1,P), sigma = 295)
initials1<- list(initials, initials, initials)

no_dist_seeds_femi_lo_hi <- stan(file = "FEMI_Rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","fg"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_femi_lo_hi, pars="lambda")
pairs(no_dist_seeds_femi_lo_hi)

### Save posterior distributions to file
save(no_dist_seeds_femi_lo_hi, file = "no_dist_seeds_femi_lo_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_femi_lo_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
femi_lo_hi <- rstan::extract(no_dist_seeds_femi_lo_hi)
acf(femi_lo_hi$lambda)

######################################################

# low intermediate initials
initials <- list(lambda=24, alpha_pler=0.07, alpha_brho=0.14, alpha_lapl=0.01,
                 alpha_femi=0.16, epsilon=rep(1,P), sigma = 181)
initials1<- list(initials, initials, initials)

no_dist_seeds_femi_lo_int <- stan(file = "FEMI_Constrained_rplot_four_species_BH_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","fg"),
                                  iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999999999, max_treedepth =50),
                                  init = initials1)

traceplot(no_dist_seeds_femi_lo_int, pars="lambda")
pairs(no_dist_seeds_femi_lo_int)

### Save posterior distributions to file
save(no_dist_seeds_femi_lo_int, file = "no_dist_seeds_femi_lo_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_femi_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
femi_lo_int <- rstan::extract(no_dist_seeds_femi_lo_int)
acf(femi_lo_int$lambda)

######################################################

# low low initials
initials <- list(lambda=2, alpha_pler=0.008, alpha_brho=0.005, alpha_lapl=-0.0003,
                 alpha_femi=0.008, epsilon=rep(1,P), sigma = 274)
initials1<- list(initials, initials, initials)

no_dist_seeds_femi_lo_lo <- stan(file = "FEMI_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","fg"),
                                 iter =5000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_femi_lo_lo, pars="lambda")
pairs(no_dist_seeds_femi_lo_lo)

### Save posterior distributions to file
save(no_dist_seeds_femi_lo_lo, file = "no_dist_seeds_femi_lo_lo.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_femi_lo_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
femi_lo_lo <- rstan::extract(no_dist_seeds_femi_lo_lo)
acf(femi_lo_lo$lambda)