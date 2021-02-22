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
intra <- brho

Plot <- dat$block
# high low initials
initials <- list(lambda=.27, alpha_pler=0.02, alpha_brho=0.003, alpha_lapl=-0.001,
                 alpha_femi=0.16, epsilon=rep(1,P), sigma = 135)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
no_dist_seeds_brho_lo_lo <- stan(file = "Constrained_rplot_four_species_BH_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot"),
                                 iter = 14000, chains = 3, thin = 4, control = list(adapt_delta = 0.99999999999999, max_treedepth =50),
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

