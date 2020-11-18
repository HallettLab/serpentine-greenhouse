## Loads rstan and sets the number of useable cores based on your computer
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read in data
data <- read.csv(paste(datpath, "/model_dat2.csv", sep = "")) %>%
  select(-X)
data[is.na(data)] <- 0

## Set initial parameter estimates, once for each chain to run
initials <- list(lambda=exp(10), alpha_pler=exp(0.03), alpha_brho=exp(0.03), alpha_lapl=exp(0.03),
                 alpha_femi=exp(0.03))
initials1<- list(initials, initials, initials)

## Subset data for competitor and treatment of interest
dat <- subset(data, species == "BRHO")
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

### Set intra-specific species
intra <- brho

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
no_dist_seeds_brho_lo_lo <- stan(file = "Four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "pler", "brho",
                                                                               "lapl", "femi"),
                                 iter = 40000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 10),
                                 init = initials1)

### Save posterior distributions to file
save(no_dist_seeds_brho_hi_hi, file = "brho_hi_hi_posteriors.rdata")

## Look at resulting estimated parameter distributions
stan_dens(no_dist_seeds_brho_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
brho_hi_hi <- rstan::extract(no_dist_seeds_brho_hi_hi)

