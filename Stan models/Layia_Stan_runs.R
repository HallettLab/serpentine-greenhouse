## Loads rstan and sets the number of useable cores based on your computer
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read in data
data <- read.csv(paste(datpath, "/model_dat2.csv", sep = "")) %>%
  select(-X)


#### LAPL Hi Hi ####
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
full_seeds_lapl_hi_hi <- stan(file = "LAPL_hi_hi_four_species_BH_model.stan", 
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

#### LAPL Hi Int ####
## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "LAPL")
dat <- subset(dat, waterN_treatment == "hi_int") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 1)
dat <- subset(dat, PLER_seeds_in == 0)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))

### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- lapl

Plot <- dat$block

lg <- 0.32

######################################################

# high high initials
initials <- list(lambda=mean(Fecundity),  epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
fecundity_lapl_hi_int <- stan(file = "LAPL_basic_BHmodel_v2.stan", 
                                 data = c("N", "Fecundity", "P", "Plot"),
                                 iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                                 init = initials1)

######################################################

## Now let's run our full Beverton Holt model with the competitors

dat_full <- subset(data, species == "LAPL")
dat_full <- subset(dat_full, waterN_treatment == "hi_int") %>%
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
initials <- list(lambda=10.77, alpha_pler=0.2, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 0
full_seeds_lapl_hi_int <- stan(file = "LAPL_hi_int_four_species_BH_model.stan", 
                              data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","lg"),
                              iter = 8000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth =50),
                              init = initials1)

traceplot(full_seeds_lapl_hi_int, pars="lambda")
#pairs(full_seeds_lapl_hi_hi)

### Save posterior distributions to file
save(full_seeds_lapl_hi_int, file = "full_seeds_lapl_hi_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(full_seeds_lapl_hi_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_hi_int <- rstan::extract(full_seeds_lapl_hi_int)
acf(lapl_hi_int$lambda)

#### LAPL Dry Int ####
## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "LAPL")
dat <- subset(dat, waterN_treatment == "lo_int") %>%
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

# lo int initials
initials <- list(lambda=mean(Fecundity),  epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
fecundity_lapl_lo_int <- stan(file = "LAPL_basic_BHmodel_v2.stan", 
                                 data = c("N", "Fecundity", "P", "Plot"),
                                 iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                                 init = initials1)

traceplot(no_dist_seeds_lapl_hi_hi, pars="lambda")
pairs(no_dist_seeds_lapl_hi_hi)

######################################################

## Now let's run our full Beverton Holt model with the competitors

dat_full <- subset(data, species == "LAPL")
dat_full <- subset(dat_full, waterN_treatment == "lo_int") %>%
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
full_seeds_lapl_lo_int <- stan(file = "LAPL_lo_int_four_species_BH_model.stan", 
                              data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","lg"),
                              iter = 8000, chains = 3, thin = 4, control = list(adapt_delta = 0.95, max_treedepth =50),
                              init = initials1)

traceplot(full_seeds_lapl_lo_int, pars="lambda")
pairs(full_seeds_lapl_lo_int)

### Save posterior distributions to file
save(full_seeds_lapl_lo_int, file = "full_seeds_lapl_lo_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(full_seeds_lapl_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_hi_hi <- rstan::extract(full_seeds_lapl_lo_int)
acf(lapl_hi_hi$lambda)

#### LAPL Dry High ####
## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "LAPL")
dat <- subset(dat, waterN_treatment == "lo_hi") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 1)
dat <- subset(dat, PLER_seeds_in == 0)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))

### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- lapl

Plot <- dat$block

lg <- 0.32

######################################################

# lo high initials
initials <- list(lambda=mean(Fecundity),  epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
fecundity_lapl_lo_hi <- stan(file = "LAPL_basic_BHmodel_v2.stan", 
                              data = c("N", "Fecundity", "P", "Plot"),
                              iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                              init = initials1)

######################################################

## Now let's run our full Beverton Holt model with the competitors

dat_full <- subset(data, species == "LAPL")
dat_full <- subset(dat_full, waterN_treatment == "lo_hi") %>%
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
initials <- list(lambda=39.9, alpha_pler=0.2, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit model
### Specify the stan model, the data to send it by name, iterations, chains (based on number of cores), 
### thinning constant (2 or 3 is usually fine), 
# NOTE: NUMBER OF DIVERGENT TRANSITIONS IS 4; ALL RHAT VALUES ARE 1.00
full_seeds_lapl_lo_hi <- stan(file = "LAPL_lo_hi_four_species_BH_model.stan", 
                               data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","lg"),
                               iter = 8000, chains = 3, thin = 4, control = list(adapt_delta = 0.95, max_treedepth =50),
                               init = initials1)

traceplot(full_seeds_lapl_lo_int, pars="lambda")
pairs(full_seeds_lapl_lo_int)

### Save posterior distributions to file
save(full_seeds_lapl_lo_int, file = "full_seeds_lapl_lo_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(full_seeds_lapl_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl", "alpha_femi"))

## Extract all parameter estimates
lapl_hi_hi <- rstan::extract(full_seeds_lapl_lo_int)
acf(lapl_hi_hi$lambda)

