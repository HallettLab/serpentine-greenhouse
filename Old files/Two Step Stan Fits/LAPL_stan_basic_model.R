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
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "LAPL")
dat <- subset(dat, waterN_treatment == "hi_hi") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 1)
dat <- subset(dat, PLER_seeds_in == 0)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

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
intra <- lapl

Plot <- dat$block

lg <- 0.32

######################################################

# high high initials
initials <- list(lambda=77.75,  epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_lapl_hi_hi <- stan(file = "LAPL_basic_BHmodel_v2.stan", 
                                 data = c("N", "Fecundity", "P", "Plot"),
                                 iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_lapl_hi_hi, pars="lambda")
pairs(no_dist_seeds_lapl_hi_hi)

### Save posterior distributions to file
save(no_dist_seeds_lapl_hi_hi, file = "no_dist_seeds_lapl_hi_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_lapl_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_hi_hi <- rstan::extract(no_dist_seeds_lapl_hi_hi)
acf(lapl_hi_hi$lambda)

######################################################

## Now let's run our full Beverton Holt model with the competitors

dat_full <- subset(data, species == "LAPL")
dat_full <- subset(dat_full, waterN_treatment == "hi_hi") %>%
  na.omit()

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat_full$seeds_out))

### Set population context for each species as seeds in
pler <- as.integer(dat_full$PLER_seeds_in)
brho <- as.integer(dat_full$BRHO_seeds_in)
lapl <- as.integer(dat_full$LAPL_seeds_in)
femi <- as.integer(dat_full$FEMI_seeds_in)

### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

intra <- lapl

Plot <- dat_full$block

lg <- 0.32

# high intermediate initials
initials <- list(lambda=77, alpha_pler=0.5, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
full_seeds_lapl_hi_hi <- stan(file = "LAPL_Constrained_rplot_four_species_BH_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","lg"),
                                  iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                                  init = initials1)

traceplot(full_seeds_lapl_hi_hi, pars="lambda")
pairs(full_seeds_lapl_hi_hi)

### Save posterior distributions to file
save(full_seeds_lapl_hi_hi, file = "full_seeds_lapl_hi_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(full_seeds_lapl_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_hi_hi <- rstan::extract(full_seeds_lapl_hi_hi)
acf(lapl_hi_hi$lambda)

######################################################

# high low initials
initials <- list(lambda=0.23, alpha_pler=3, alpha_brho=4, alpha_lapl=4,
                 alpha_femi=3, epsilon=rep(1,P), sigma = 197)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_lapl_hi_lo <- stan(file = "LAPL_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","lg"),
                                 iter = 10000, chains = 3, thin = 5, control = list(adapt_delta = 0.9999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_lapl_hi_lo, pars="lambda")
pairs(no_dist_seeds_lapl_hi_lo)

### Save posterior distributions to file
save(no_dist_seeds_lapl_hi_lo, file = "no_dist_seeds_lapl_hi_lo")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_lapl_hi_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_hi_lo <- rstan::extract(no_dist_seeds_lapl_hi_lo)
acf(lapl_hi_lo$lambda)

######################################################

# low high initials
initials <- list(lambda=254, alpha_pler=0.23, alpha_brho=2.3, alpha_lapl=0.8,
                 alpha_femi=1, epsilon=rep(1,P), sigma = 17)
initials1<- list(initials, initials, initials)

no_dist_seeds_lapl_lo_hi <- stan(file = "LAPL_Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","lg"),
                                 iter = 5000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_lapl_lo_hi, pars="lambda")
pairs(no_dist_seeds_lapl_lo_hi)

### Save posterior distributions to file
save(no_dist_seeds_lapl_lo_hi, file = "no_dist_seeds_lapl_lo_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_lapl_lo_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_lo_hi <- rstan::extract(no_dist_seeds_lapl_lo_hi)
acf(lapl_lo_hi$lambda)

######################################################

# low intermediate initials
initials <- list(lambda=227, alpha_pler=7.5, alpha_brho=7.6, alpha_lapl=4.8,
                 alpha_femi=6.5, epsilon=rep(1,P), sigma = 1.8)
initials1<- list(initials, initials, initials)

no_dist_seeds_lapl_lo_int <- stan(file = "LAPL_Constrained_rplot_four_species_BH_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","lg"),
                                  iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999999999999999, max_treedepth =50),
                                  init = initials1)

traceplot(no_dist_seeds_lapl_lo_int, pars="lambda")
pairs(no_dist_seeds_lapl_lo_int)

### Save posterior distributions to file
save(no_dist_seeds_lapl_lo_int, file = "no_dist_seeds_lapl_lo_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_lapl_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_lo_int <- rstan::extract(no_dist_seeds_lapl_lo_int)
acf(lapl_lo_int$lambda)

######################################################

# low low initials
initials <- list(lambda=1, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 alpha_femi=0.05, epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

no_dist_seeds_lapl_lo_lo <- stan(file = "Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot"),
                                 iter =1000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999999999999999, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_lapl_lo_lo, pars="lambda")
pairs(no_dist_seeds_lapl_lo_lo)

### Save posterior distributions to file
save(no_dist_seeds_lapl_lo_lo, file = "no_dist_seeds_lapl_lo_lo")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_lapl_lo_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_lo_lo <- rstan::extract(no_dist_seeds_lapl_lo_lo)
acf(lapl_lo_lo$lambda)