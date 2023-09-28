## Loads rstan and sets the number of useable cores based on your computer
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read in data
data <- read.csv(paste(datpath, "model_dat2_3spp.csv", sep = "")) %>%
  select(-X) 

## Subset data for competitor and treatment of interest
dat <- subset(data, species == "LAPL")
dat <- subset(dat, waterN_treatment == "hi_hi") %>%
  na.omit()

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))

### Set population context for each species as seeds in
pler <- as.integer(dat$PLER_seeds_in)
brho <- as.integer(dat$BRHO_seeds_in)
lapl <- as.integer(dat$LAPL_seeds_in)


### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- lapl
Plot <- dat$block

lg <- 0.32
bg <- 0.98
pg <- 0.92

######################################################
# HIGH HIGH
# high high dirty run
initials <- list(lambda=78, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

lapl_hi_hi <- stan(file = "lapl_constrained_bevertonholt_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl","P", "Plot","lg","bg","pg"),
                                 iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                                 init = initials1)

# high high model fit
initials <- list(lambda=602, alpha_pler=0.43, alpha_brho=10.49, alpha_lapl=4.34,
                 epsilon=rep(1,P), sigma = 5)
initials1<- list(initials, initials, initials)
lapl_hi_hi <- stan(file = "lapl_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl","P", "Plot","lg","bg","pg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 10),
                   init = initials1)

traceplot(lapl_hi_hi, pars="lambda")
pairs(lapl_hi_hi)

### Save posterior distributions to file
save(lapl_hi_hi, file = "lapl_hi_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(lapl_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
lapl_hi_hi <- rstan::extract(lapl_hi_hi)
acf(lapl_hi_hi$lambda)

######################################################
# HIGH INTERMEDIATE
## Subset data for competitor and treatment of interest
dat <- subset(data, species == "LAPL")
dat <- subset(dat, waterN_treatment == "hi_int") %>%
  na.omit()

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))

### Set population context for each species as seeds in
pler <- as.integer(dat$PLER_seeds_in)
brho <- as.integer(dat$BRHO_seeds_in)
lapl <- as.integer(dat$LAPL_seeds_in)


### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- lapl
Plot <- dat$block

# high intermediate dirty run
initials <- list(lambda=11, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma =10)
initials1<- list(initials, initials, initials)

lapl_hi_int <- stan(file = "lapl_constrained_bevertonholt_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","lg","bg","pg"),
                                  iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =10),
                                  init = initials1)

# high intermediate model fit
initials <- list(lambda=117, alpha_pler=0.34, alpha_brho=6.61, alpha_lapl=9.86,
                 epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

lapl_hi_int <- stan(file = "lapl_bevertonholt_model.stan", 
                    data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","lg","bg","pg"),
                    iter = 20000, chains = 3, thin = 10, control = list(adapt_delta = 0.99999999999, max_treedepth =50),
                    init = initials1)

traceplot(lapl_hi_int, pars="lambda")
pairs(lapl_hi_int)

### Save posterior distributions to file
save(lapl_hi_int, file = "lapl_hi_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(lapl_hi_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
lapl_hi_int <- rstan::extract(lapl_hi_int)
acf(lapl_hi_int$lambda)

######################################################
# HIGH LOW
# NA - no seeds produced

######################################################
# LOW HIGH
## Subset data for competitor and treatment of interest
dat <- subset(data, species == "LAPL")
dat <- subset(dat, waterN_treatment == "lo_hi") %>%
  na.omit()

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))

### Set population context for each species as seeds in
pler <- as.integer(dat$PLER_seeds_in)
brho <- as.integer(dat$BRHO_seeds_in)
lapl <- as.integer(dat$LAPL_seeds_in)


### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- lapl
Plot <- dat$block

# low high dirty run
initials <- list(lambda=40, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

lapl_lo_hi <- stan(file = "lapl_constrained_bevertonholt_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl","P", "Plot","lg","bg","pg"),
                                 iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =10),
                                 init = initials1)

# low high model fit
initials <- list(lambda=245, alpha_pler=0.24, alpha_brho=2.33, alpha_lapl=2.52,
                 epsilon=rep(1,P), sigma = 18)
initials1<- list(initials, initials, initials)

lapl_lo_hi <- stan(file = "lapl_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl","P", "Plot","lg","bg","pg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth =10),
                   init = initials1)

traceplot(lapl_lo_hi, pars="lambda")
pairs(lapl_lo_hi)

### Save posterior distributions to file
save(lapl_lo_hi, file = "lapl_lo_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(lapl_lo_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
lapl_lo_hi <- rstan::extract(lapl_lo_hi)
acf(lapl_lo_hi$lambda)

######################################################
# LOW INTERMEDIATE
## Subset data for competitor and treatment of interest
dat <- subset(data, species == "LAPL")
dat <- subset(dat, waterN_treatment == "lo_int") %>%
  na.omit()

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))

### Set population context for each species as seeds in
pler <- as.integer(dat$PLER_seeds_in)
brho <- as.integer(dat$BRHO_seeds_in)
lapl <- as.integer(dat$LAPL_seeds_in)


### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- lapl
Plot <- dat$block

# low intermediate dirty run
initials <- list(lambda=17, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

lapl_lo_int <- stan(file = "lapl_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl","P", "Plot","lg","bg","pg"),
                   iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =10),
                   init = initials1)

# low intermediate model fit
initials <- list(lambda=245, alpha_pler=0.24, alpha_brho=2.33, alpha_lapl=2.52,
                 epsilon=rep(1,P), sigma = 18)
initials1<- list(initials, initials, initials)

lapl_lo_int <- stan(file = "lapl_bevertonholt_model.stan", 
                    data = c("N", "Fecundity", "intra", "pler", "brho", "lapl","P", "Plot","lg","bg","pg"),
                    iter = 20000, chains = 3, thin = 10, control = list(adapt_delta = 0.9999999999, max_treedepth = 50),
                    init = initials1)

traceplot(no_dist_seeds_lapl_lo_int, pars="lambda")
pairs(no_dist_seeds_lapl_lo_int)

### Save posterior distributions to file
save(lapl_lo_int, file = "lapl_lo_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(lapl_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
lapl_lo_int <- rstan::extract(lapl_lo_int)
acf(lapl_lo_int$lambda)

######################################################
# Low LOW
# NA - no seeds produced