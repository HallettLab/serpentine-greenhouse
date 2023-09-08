## Loads rstan and sets the number of useable cores based on your computer
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read in data
data <- read.csv(paste(datpath, "/model_dat2_3spp.csv", sep = "")) %>%
  select(-X) %>%
  mutate(BRHO_seeds_in = ifelse(BRHO_seeds_in %in% 1, 1/0.98, BRHO_seeds_in)) %>%
  mutate(PLER_seeds_in = ifelse(PLER_seeds_in %in% 1, 1/0.92, PLER_seeds_in)) %>%
  mutate(LAPL_seeds_in = ifelse(LAPL_seeds_in %in% 1, 1/0.32, LAPL_seeds_in))

## Subset data for competitor and treatment of interest
dat <- subset(data, species == "BRHO")
dat <- subset(dat, waterN_treatment == "hi_hi") %>%
  na.omit()

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))

### Set population context for each species as seeds in
pler <- dat$PLER_seeds_in
brho <- dat$BRHO_seeds_in
lapl <- dat$LAPL_seeds_in

### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- brho

Plot <- dat$block

bg <- 0.98
pg <- 0.92
lg <- 0.32

######################################################
# HIGH HIGH
# high high dirty run
initials <- list(lambda=113, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1 <- list(initials, initials, initials)

brho_hi_hi <- stan(file = "brho_bevertonholt_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                                 iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =10),
                   init = initials1)

# high high model fit
initials <- list(lambda=240, alpha_pler=0.24, alpha_brho=1.11, alpha_lapl=0.74,
                 epsilon=rep(1,P), sigma = 296)
initials1 <- list(initials, initials, initials)

brho_hi_hi <- stan(file = "brho_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =10),
                   init = initials1)

traceplot(brho_hi_hi, pars="lambda")
pairs(brho_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

### Save posterior distributions to file
save(brho_hi_hi, file = "brho_hi_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(brho_hi_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
brho_hi_hi <- rstan::extract(brho_hi_hi)
acf(brho_hi_hi$lambda)

######################################################
# HIGH INTERMEDIATE
## Subset data for competitor and treatment of interest
dat <- subset(data, species == "BRHO")
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
intra <- brho

Plot <- dat$block

# high intermediate dirty run
initials <- list(lambda=28, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

brho_hi_int <- stan(file = "brho_bevertonholt_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                                 iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =10),
                                 init = initials1)

# high intermediate model fit
initials <- list(lambda=51, alpha_pler=0.17, alpha_brho=0.81, alpha_lapl=0.56,
                 epsilon=rep(1,P), sigma = 163)
initials1<- list(initials, initials, initials)

brho_hi_int <- stan(file = "brho_bevertonholt_model.stan", 
                    data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                    iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =10),
                    init = initials1)

traceplot(brho_hi_lo, pars="lambda")
pairs(brho_hi_int)

### Save posterior distributions to file
save(brho_hi_int, file = "brho_hi_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(brho_hi_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
brho_hi_int <- rstan::extract(brho_hi_int)
acf(brho_hi_int$lambda)

######################################################
# HIGH LOW - LAUREN S HELP?
## Subset data for competitor and treatment of interest
dat <- subset(data, species == "BRHO")
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
intra <- brho

Plot <- dat$block

# high low dirty run
initials <- list(lambda=7.5, alpha_pler=0.1, alpha_brho=0.3, alpha_lapl=0.1,
                 epsilon=rep(1,P), sigma = 220)
initials <- list(lambda=7, alpha_pler=0.1, alpha_brho=0.3, alpha_lapl=0.1)
initials1<- list(initials, initials, initials)

brho_hi_lo <- stan(file = "brho_hi_lo_bevertonholt_model.stan", 
                    data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                   iter = 20000, chains = 3, thin = 12, control = list(adapt_delta = 0.9999999999999999, max_treedepth =100, stepsize=0.1),
                    init = initials1, cores=3)

brho_hi_lo <- stan(file = "brho_hi_lo_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.999, max_treedepth =20),
                   init = initials1, cores=3)

# high low model fit
#initials <- list(lambda=7.83, alpha_pler=0.11, alpha_brho=0.38, alpha_lapl=0.15,
#                 epsilon=rep(1,P), sigma = 242)
#initials1<- list(initials, initials, initials)

#brho_hi_lo <- stan(file = "brho_bevertonholt_model.stan", 
#                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
#                   iter = 16000, warmup=10000, chains = 3, thin = 3, control = list(adapt_delta = 0.9999, max_treedepth =100),
#                                 init = initials1)

traceplot(brho_hi_lo, pars="lambda")
pairs(brho_hi_lo)

### Save posterior distributions to file
save(brho_hi_lo, file = "brho_hi_lo.rdata")

## Look at resulting estimated parameter distributions
stan_dens(brho_hi_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
brho_hi_lo <- rstan::extract(brho_hi_lo)
acf(brho_hi_lo$lambda)

######################################################
# LOW HIGH - LAUREN S HELP?
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

### Number of observations 
N <- as.integer(length(Fecundity))
P <- as.integer(length(unique(dat$block)))

### Set intra-specific species
intra <- brho
Plot <- dat$block

# lo hi dirty run
initials <- list(lambda=85, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

brho_lo_hi <- stan(file = "brho_lo_hi_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                   iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.999, max_treedepth = 10),
                   init = initials1)

#low high model fit
initials <- list(lambda=273, alpha_pler=0.02, alpha_brho=0.16, alpha_lapl=0.07,
                 epsilon=rep(1,P), sigma = 147)
initials1<- list(initials, initials, initials)

brho_lo_hi <- stan(file = "brho_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                   iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 20),
                   init = initials1)


traceplot(brho_lo_hi, pars="lambda")
pairs(brho_lo_hi)

### Save posterior distributions to file
save(brho_lo_hi, file = "brho_lo_hi.rdata")

## Look at resulting estimated parameter distributions
stan_dens(brho_lo_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

pairs(brho_lo_hi, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
brho_lo_hi <- rstan::extract(brho_lo_hi)
cor(brho_lo_hi$lambda, brho_lo_hi$alpha_pler) #-.03
cor(brho_lo_hi$lambda, brho_lo_hi$alpha_lapl) #.003
cor(brho_lo_hi$lambda, brho_lo_hi$alpha_brho) #-0.012

acf(brho_lo_hi$lambda)

######################################################
# LOW INTERMEDIATE
dat <- subset(data, species == "BRHO")
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
intra <- brho
Plot <- dat$block

# low intermediate dirty run
initials <- list(lambda=14, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl=0.05,
                 epsilon=rep(1,P), sigma = 18)
initials1<- list(initials, initials, initials)

brho_lo_int <- stan(file = "brho_constrained_bevertonholt_model.stan", 
                                 data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                                 iter = 2000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                                 init = initials1)

# low intermediate model fit
initials <- list(lambda=24, alpha_pler=0.08, alpha_brho=0.54, alpha_lapl=0.13,
                 epsilon=rep(1,P), sigma = 33)
initials1<- list(initials, initials, initials)

brho_lo_int <- stan(file = "brho_constrained_bevertonholt_model.stan", 
                    data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "P", "Plot","bg","pg","lg"),
                    iter = 12000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth = 10),
                    init = initials1)

traceplot(brho_lo_int, pars="lambda")
pairs(brho_lo_int)

### Save posterior distributions to file
save(brho_lo_int, file = "brho_lo_int.rdata")

## Look at resulting estimated parameter distributions
stan_dens(brho_lo_int, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
brho_lo_int <- rstan::extract(brho_lo_int)
acf(brho_lo_int$lambda)

######################################################
# LOW LOW - LAUREN S HELP?
dat <- subset(data, species == "BRHO")
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
intra <- brho
Plot <- dat$block

# low low dirty run
initials <- list(lambda=2.25, alpha_pler=0.05, alpha_brho=0.05, alpha_lapl= 0.05,
                 epsilon=rep(1,P), sigma = 10)
initials1<- list(initials, initials, initials)

brho_lo_lo <- stan(file = "brho2_constrained_bevertonholt_model.stan", 
                                  data = c("N", "Fecundity", "intra", "pler", "brho", "lapl","P", "Plot","bg","pg","lg"),
                                  iter = 12000, chains = 3, thin = 5, control = list(adapt_delta = 0.9, max_treedepth = 10),
                                  init = initials1)

# low low model fit
initials <- list(lambda=0.24, alpha_pler=0.01, alpha_brho=0.0006, alpha_lapl= -0.005,
                 epsilon=rep(1,P), sigma = 160)
initials1<- list(initials, initials, initials)

brho_lo_lo <- stan(file = "brho2_constrained_bevertonholt_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl","P", "Plot","bg","pg","lg"),
                   iter = 12000, chains = 3, thin = 5, control = list(adapt_delta = 0.9, max_treedepth = 20),
                   init = initials1)


traceplot(brho_lo_lo, pars="lambda")
pairs(brho_lo_lo)

### Save posterior distributions to file
save(brho_lo_lo, file = "brho_lo_lo.rdata")

## Look at resulting estimated parameter distributions
stan_dens(brho_lo_lo, pars = c("lambda", "alpha_pler", "alpha_brho", "alpha_lapl"))

## Extract all parameter estimates
brho_lo_lo <- rstan::extract(brho_lo_lo)
acf(brho_lo_lo$lambda)