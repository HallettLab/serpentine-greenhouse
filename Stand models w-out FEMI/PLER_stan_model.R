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
dat <- subset(data, species == "PLER")
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
intra <- pler

Plot <- dat$block

pg <- 0.92

######################################################
# high high initials
initials <- list(lambda=15.7, alpha_pler=0.10, alpha_brho=0.14, alpha_lapl=0.03,
                 alpha_femi=0.10, epsilon=rep(1,P), sigma = 12.6)
initials1<- list(initials, initials, initials)

## Fit high high high model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_pler_hi_hi <- stan(file = "PLER_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_pler_hi_hi, pars="lambda")
pairs(no_dist_seeds_pler_hi_hi)

### Save posterior distributions to file
save(no_dist_seeds_pler_hi_hi, file = "no_dist_seeds_pler_hi_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_pler_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
pler_hi_hi <- rstan::extract(no_dist_seeds_pler_hi_hi)
acf(pler_hi_hi$lambda)

#############################################################

# High intermediate initials
initials <- list(lambda=5.5, alpha_pler=0.08, alpha_brho=0.09, alpha_lapl=0.03,
                 alpha_femi=0.02, epsilon=rep(1,P), sigma = 21.4)
initials1<- list(initials, initials, initials)

## Fit high high high model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_pler_hi_int <- stan(file = "PLER_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_pler_hi_int, pars="lambda")
pairs(no_dist_seeds_pler_hi_int)

### Save posterior distributions to file
save(no_dist_seeds_pler_hi_int, file = "no_dist_seeds_pler_hi_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_pler_hi_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
pler_hi_int <- rstan::extract(no_dist_seeds_pler_hi_int)
acf(pler_hi_int$lambda)

######################################################
# High low initials
initials <- list(lambda=1.8, alpha_pler=0.01, alpha_brho=0.005, alpha_lapl=0.007,
                 alpha_femi=0.005, epsilon=rep(1,P), sigma = 4.8)
initials1<- list(initials, initials, initials)

no_dist_seeds_pler_hi_lo <- stan(file = "PLER_Constrained_rplot_four_species_BH_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg"),
                                  iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999999999999999, max_treedepth =50),
                                  init = initials1)

traceplot(no_dist_seeds_pler_hi_lo, pars="lambda")
pairs(no_dist_seeds_pler_hi_lo)

### Save posterior distributions to file
save(no_dist_seeds_pler_hi_lo, file = "no_dist_seeds_pler_hi_lo.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_pler_hi_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
pler_hi_lo <- rstan::extract(no_dist_seeds_pler_hi_lo)
acf(pler_hi_lo$lambda)

########################################################
# Low high initials
initials <- list(lambda=29.8, alpha_pler=0.21, alpha_brho=0.27, alpha_lapl=0.08,
                 alpha_femi=0.09, epsilon=rep(1,P), sigma = 6.3)
initials1<- list(initials, initials, initials)

no_dist_seeds_pler_lo_hi <- stan(file = "PLER_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.999999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_pler_lo_hi, pars="lambda")
pairs(no_dist_seeds_pler_lo_hi)

### Save posterior distributions to file
save(no_dist_seeds_pler_lo_hi, file = "no_dist_seeds_pler_lo_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_pler_lo_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
pler_lo_hi <- rstan::extract(no_dist_seeds_pler_lo_hi)
acf(pler_lo_hi$lambda)

######################################################

# Low intermediate initials
initials <- list(lambda=8.4, alpha_pler=0.09, alpha_brho=0.06, alpha_lapl=0.02,
                 alpha_femi=0.04, epsilon=rep(1,P), sigma = 3.2)
initials1<- list(initials, initials, initials)

no_dist_seeds_pler_lo_int <- stan(file = "PLER_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg"),
                                 iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_pler_lo_int, pars="lambda")
pairs(no_dist_seeds_pler_lo_int)

### Save posterior distributions to file
save(no_dist_seeds_pler_lo_int, file = "no_dist_seeds_pler_lo_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_pler_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
pler_lo_int <- rstan::extract(no_dist_seeds_pler_lo_int)
acf(pler_lo_int$lambda)

########################################################

# Low low initials
initials <- list(lambda=2.04, alpha_pler=0.029, alpha_brho=0.03, alpha_lapl=0.004,
                 alpha_femi=0.014, epsilon=rep(1,P), sigma = 4.7)
initials1<- list(initials, initials, initials)

no_dist_seeds_pler_lo_lo <- stan(file = "PLER_Constrained_rplot_four_species_BH_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg"),
                                  iter = 10000, chains = 3, thin = 3, control = list(adapt_delta = 0.999999999999999, max_treedepth =50),
                                  init = initials1)

traceplot(no_dist_seeds_pler_lo_lo, pars="lambda")
pairs(no_dist_seeds_pler_lo_lo)

### Save posterior distributions to file
save(no_dist_seeds_pler_lo_lo, file = "no_dist_seeds_pler_lo_lo.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_pler_lo_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
pler_lo_lo <- rstan::extract(no_dist_seeds_pler_lo_lo)
acf(pler_lo_lo$lambda)