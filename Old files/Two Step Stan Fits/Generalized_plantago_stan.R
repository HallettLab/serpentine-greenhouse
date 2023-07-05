## Loads rstan and sets the number of useable cores based on your computer
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read in data
data <- read.csv(paste(datpath, "/model_dat2.csv", sep = "")) %>%
  select(-X)

#### PLER Hi Hi -------------------------------------------------------------
#### a. Fecundity Fit ####

## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "hi_hi") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 0)
dat <- subset(dat, PLER_seeds_in == 1)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))
N <- as.integer(length(Fecundity))

# set initial starting values
initials <- list(lambda=mean(Fecundity))
initials1<- list(initials, initials, initials)

## Fit basic fecundity model
fecundity_pler_hi_hi <- stan(file = "PLER_fecundity_model.stan", 
                             data = c("N", "Fecundity"),
                             iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                             init = initials1)

#### b. Full Fit ####
## Now run full Beverton Holt model with the competitors
dat_full <- subset(data, species == "PLER")
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
intra <- pler
Plot <- dat_full$block
pg <- 0.92

# Set fecundity posterior
fecund_fit <- summary(fecundity_pler_hi_hi)$summary
fecund_mean <- fecund_fit["lambda", "mean"]
fecund_sd <- fecund_fit["lambda", "sd"]

initials <- list(lambda=fecund_mean, alpha_pler=0.2, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit full model
pler_hi_hi <- stan(file = "PLER_four_species_BH_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg", "fecund_mean", "fecund_sd"),
                   iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                   init = initials1)

### Save posterior distributions to file
save(pler_hi_hi, file = "pler_hi_hi.rdata")

#### PLER Hi Int -------------------------------------------------------------
#### a. Fecundity Fit ####

## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "hi_int") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 0)
dat <- subset(dat, PLER_seeds_in == 1)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))
N <- as.integer(length(Fecundity))

# set initial starting values
initials <- list(lambda=mean(Fecundity))
initials1<- list(initials, initials, initials)

## Fit basic fecundity model
fecundity_pler_hi_int <- stan(file = "PLER_fecundity_model.stan", 
                             data = c("N", "Fecundity"),
                             iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                             init = initials1)

#### b. Full Fit ####
## Now run full Beverton Holt model with the competitors
dat_full <- subset(data, species == "PLER")
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
intra <- pler
Plot <- dat_full$block
pg <- 0.92

# Set fecundity posterior
fecund_fit <- summary(fecundity_pler_hi_int)$summary
fecund_mean <- fecund_fit["lambda", "mean"]
fecund_sd <- fecund_fit["lambda", "sd"]

initials <- list(lambda=fecund_mean, alpha_pler=0.2, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit full model
pler_hi_int <- stan(file = "PLER_four_species_BH_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg", "fecund_mean", "fecund_sd"),
                   iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth =50),
                   init = initials1)

### Save posterior distributions to file
save(pler_hi_int, file = "pler_hi_int.rdata")

#### PLER Hi LO -------------------------------------------------------------
#### a. Fecundity Fit ####

## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "hi_lo") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 0)
dat <- subset(dat, PLER_seeds_in == 1)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))
N <- as.integer(length(Fecundity))

# set initial starting values
initials <- list(lambda=mean(Fecundity))
initials1<- list(initials, initials, initials)

## Fit basic fecundity model
fecundity_pler_hi_lo <- stan(file = "PLER_fecundity_model.stan", 
                             data = c("N", "Fecundity"),
                             iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                             init = initials1)

#### b. Full Fit ####
## Now run full Beverton Holt model with the competitors
dat_full <- subset(data, species == "PLER")
dat_full <- subset(dat_full, waterN_treatment == "hi_lo") %>%
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
intra <- pler
Plot <- dat_full$block
pg <- 0.92

# Set fecundity posterior
fecund_fit <- summary(fecundity_pler_hi_lo)$summary
fecund_mean <- fecund_fit["lambda", "mean"]
fecund_sd <- fecund_fit["lambda", "sd"]

initials <- list(lambda=fecund_mean, alpha_pler=0.2, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit full model
pler_hi_lo <- stan(file = "PLER_four_species_BH_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg", "fecund_mean", "fecund_sd"),
                   iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth =50),
                   init = initials1)

### Save posterior distributions to file
save(pler_hi_lo, file = "pler_hi_lo.rdata")

#### PLER Lo Hi -------------------------------------------------------------
#### a. Fecundity Fit ####

## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "lo_hi") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 0)
dat <- subset(dat, PLER_seeds_in == 1)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))
N <- as.integer(length(Fecundity))

# set initial starting values
initials <- list(lambda=mean(Fecundity))
initials1<- list(initials, initials, initials)

## Fit basic fecundity model
fecundity_pler_lo_hi <- stan(file = "PLER_fecundity_model.stan", 
                             data = c("N", "Fecundity"),
                             iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                             init = initials1)

#### b. Full Fit ####
## Now run full Beverton Holt model with the competitors
dat_full <- subset(data, species == "PLER")
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
intra <- pler
Plot <- dat_full$block
pg <- 0.92

# Set fecundity posterior
fecund_fit <- summary(fecundity_pler_lo_hi)$summary
fecund_mean <- fecund_fit["lambda", "mean"]
fecund_sd <- fecund_fit["lambda", "sd"]

initials <- list(lambda=fecund_mean, alpha_pler=0.2, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit full model
pler_lo_hi <- stan(file = "PLER_four_species_BH_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg", "fecund_mean", "fecund_sd"),
                   iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                   init = initials1)

### Save posterior distributions to file
save(pler_lo_hi, file = "pler_lo_hi.rdata")

#### PLER Lo Int -------------------------------------------------------------
#### a. Fecundity Fit ####

## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "lo_int") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 0)
dat <- subset(dat, PLER_seeds_in == 1)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))
N <- as.integer(length(Fecundity))

# set initial starting values
initials <- list(lambda=mean(Fecundity))
initials1<- list(initials, initials, initials)

## Fit basic fecundity model
fecundity_pler_lo_int <- stan(file = "PLER_fecundity_model.stan", 
                             data = c("N", "Fecundity"),
                             iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                             init = initials1)

#### b. Full Fit ####
## Now run full Beverton Holt model with the competitors
dat_full <- subset(data, species == "PLER")
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
intra <- pler
Plot <- dat_full$block
pg <- 0.92

# Set fecundity posterior
fecund_fit <- summary(fecundity_pler_lo_int)$summary
fecund_mean <- fecund_fit["lambda", "mean"]
fecund_sd <- fecund_fit["lambda", "sd"]

initials <- list(lambda=fecund_mean, alpha_pler=0.2, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit full model
pler_lo_int <- stan(file = "PLER_four_species_BH_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg", "fecund_mean", "fecund_sd"),
                   iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                   init = initials1)

### Save posterior distributions to file
save(pler_lo_int, file = "pler_lo_int.rdata")

#### PLER Lo Lo -------------------------------------------------------------
#### a. Fecundity Fit ####

## Subset data for competitor and treatment of interest
## First, we want to estimate lambda without any intra or interspecific competition
dat <- subset(data, species == "PLER")
dat <- subset(dat, waterN_treatment == "lo_lo") %>%
  na.omit()
dat <- subset(dat, LAPL_seeds_in == 0)
dat <- subset(dat, PLER_seeds_in == 1)
dat <- subset(dat, FEMI_seeds_in == 0)
dat <- subset(dat, BRHO_seeds_in == 0)

## Create model variables for our data
### Set Fecundity as the seeds out from our focal species
Fecundity <- as.integer(round(dat$seeds_out))
N <- as.integer(length(Fecundity))

# set initial starting values
initials <- list(lambda=mean(Fecundity))
initials1<- list(initials, initials, initials)

## Fit basic fecundity model
fecundity_pler_lo_lo <- stan(file = "PLER_fecundity_model.stan", 
                             data = c("N", "Fecundity"),
                             iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.9, max_treedepth =50),
                             init = initials1)

#### b. Full Fit ####
## Now run full Beverton Holt model with the competitors
dat_full <- subset(data, species == "PLER")
dat_full <- subset(dat_full, waterN_treatment == "lo_lo") %>%
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
intra <- pler
Plot <- dat_full$block
pg <- 0.92

# Set fecundity posterior
fecund_fit <- summary(fecundity_pler_lo_lo)$summary
fecund_mean <- fecund_fit["lambda", "mean"]
fecund_sd <- fecund_fit["lambda", "sd"]

initials <- list(lambda=fecund_mean, alpha_pler=0.2, alpha_brho=.2, alpha_lapl=.2,
                 alpha_femi=.2, epsilon=rep(1,P), sigma =2)
initials1<- list(initials, initials, initials)

## Fit full model
pler_lo_lo <- stan(file = "PLER_four_species_BH_model.stan", 
                   data = c("N", "Fecundity", "intra", "pler", "brho", "lapl", "femi", "P", "Plot","pg", "fecund_mean", "fecund_sd"),
                   iter = 6000, chains = 3, thin = 3, control = list(adapt_delta = 0.99, max_treedepth =50),
                   init = initials1)

### Save posterior distributions to file
save(pler_lo_lo, file = "pler_lo_lo.rdata")
