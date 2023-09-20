## Loads rstan and sets the number of useable cores based on your computer
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read in data
data <- read.csv(paste(datpath, "model_dat2_3spp.csv", sep = "")) %>%
  select(-X) %>%
  mutate(BRHO_seeds_in = ifelse(BRHO_seeds_in %in% 1, 1/0.98, BRHO_seeds_in)) %>%
  mutate(PLER_seeds_in = ifelse(PLER_seeds_in %in% 1, 1/0.92, PLER_seeds_in)) %>%
  mutate(LAPL_seeds_in = ifelse(LAPL_seeds_in %in% 1, 1/0.32, LAPL_seeds_in))

## Subset data for competitor and treatment of interest
dat <- subset(data, species == "PLER")
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
intra <- pler
Plot <- dat$block

bg <- 0.98
pg <- 0.92
lg <- 0.32

######################################################
# HIGH HIGH
# high high dirty run
initials <- list(lambda=10.5, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

pler_hi_hi <- stan(file = "pler_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                   init = initials1)

# high high model fit
initials <- list(lambda= 16, alpha_pler= 0.11, alpha_brho=0.15, alpha_lapl=0.08,
                 epsilon=rep(1,P), sigma = 13)
initials1<- list(initials, initials, initials)

pler_hi_hi <- stan(file = "pler_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                   init = initials1)


traceplot(pler_hi_hi, pars="lambda")
pairs(pler_hi_hi)

### Save posterior distributions to file
save(pler_hi_hi, file = "pler_hi_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(pler_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
pler_hi_hi <- rstan::extract(pler_hi_hi)
acf(pler_hi_hi$lambda)

#############################################################
# HIGH INTERMEDIATE
## Subset data for competitor and treatment of interest
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "hi_int")

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
intra <- pler
Plot <- dat$block

# High intermediate dirty run
initials <- list(lambda=3.5, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

pler_hi_int <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                   init = initials1)

# High intermediate model fit
initials <- list(lambda=5, alpha_pler=0.09, alpha_brho=0.09, alpha_lapl=0.10,
                 epsilon=rep(1,P), sigma = 29)
initials1<- list(initials, initials, initials)

pler_hi_int <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 10),
                   init = initials1)

traceplot(pler_hi_int, pars="lambda")
pairs(pler_hi_int)

### Save posterior distributions to file
save(pler_hi_int, file = "pler_hi_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(pler_hi_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
pler_hi_int <- rstan::extract(pler_hi_int)
acf(pler_hi_int$lambda)

######################################################
# HIGH LOW
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "hi_lo")

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
intra <- pler
Plot <- dat$block

# High low dirty run
initials <- list(lambda=1.8, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

pler_hi_lo <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                    data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                    iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                    init = initials1)

# High low model fit
initials <- list(lambda=1.6, alpha_pler=0.01, alpha_brho=0.01, alpha_lapl=0.02,
                 epsilon=rep(1,P), sigma = 4)
initials1<- list(initials, initials, initials)

pler_hi_lo <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.999, max_treedepth = 10),
                   init = initials1)


traceplot(pler_hi_lo, pars="lambda")
pairs(pler_hi_lo)

### Save posterior distributions to file
save(pler_hi_lo, file = "pler_hi_lo.rdata")

## Look at resulting estimated parameter distributions
stan_dens(pler_hi_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
pler_hi_lo <- rstan::extract(pler_hi_lo)
acf(pler_hi_lo$lambda)

########################################################
# LOW HIGH
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "lo_hi")

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
intra <- pler
Plot <- dat$block

# Low high dirty run
initials <- list(lambda= 23, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

pler_lo_hi <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                   init = initials1)

# Low high model fit
initials <- list(lambda= 32, alpha_pler=0.22, alpha_brho=0.26, alpha_lapl=0.24,
                 epsilon=rep(1,P), sigma = 6)
initials1<- list(initials, initials, initials)

pler_lo_hi <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 10),
                   init = initials1)

traceplot(pler_lo_hi, pars="lambda")
pairs(pler_lo_hi)

### Save posterior distributions to file
save(pler_lo_hi, file = "pler_lo_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(pler_lo_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
pler_lo_hi <- rstan::extract(pler_lo_hi)
acf(pler_lo_hi$lambda)

######################################################
# LOW INTERMEDIATE
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "lo_int")

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
intra <- pler
Plot <- dat$block

# Low intermediate dirty run
initials <- list(lambda= 6, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

pler_lo_int <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                   init = initials1)

# Low intermediate model fit
initials <- list(lambda= 8, alpha_pler=0.1, alpha_brho=0.06, alpha_lapl=0.06,
                 epsilon=rep(1,P), sigma = 3)
initials1<- list(initials, initials, initials)

pler_lo_int <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                    data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                    iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                    init = initials1)

traceplot(pler_lo_int, pars="lambda")
pairs(pler_lo_int)

### Save posterior distributions to file
save(pler_lo_int, file = "pler_lo_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(pler_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
pler_lo_int <- rstan::extract(pler_lo_int)
acf(pler_lo_int$lambda)

########################################################
# LOW LOW
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "lo_lo")

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
intra <- pler
Plot <- dat$block

# Low low dirty run
initials <- list(lambda= 2, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

pler_lo_lo <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                    data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                    iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                    init = initials1)

# Low low model fit
initials <- list(lambda= 2, alpha_pler=0.03, alpha_brho=0.03, alpha_lapl=0.01,
                 epsilon=rep(1,P), sigma = 5)
initials1<- list(initials, initials, initials)

pler_lo_lo <- stan(file = "pler_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","pg","bg","lg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999, max_treedepth = 30),
                   init = initials1)

traceplot(pler_lo_lo, pars="lambda")
pairs(pler_lo_lo)

### Save posterior distributions to file
save(pler_lo_lo, file = "pler_lo_lo.rdata")

## Look at resulting estimated parameter distributions
stan_dens(pler_lo_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
pler_lo_lo <- rstan::extract(pler_lo_lo)
acf(pler_lo_lo$lambda)