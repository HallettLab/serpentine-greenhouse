##########################################
####Equilibrium abundance HI H20 Hi N#####
##########################################

library(rstan)
library(here)

##survival and germination fractions
ps <- .75 # gulmon
pg <- .92 # gulmon
bs <- .013 # andrew
bg <- .98 # gulmon
ls <- .15 # rossington
lg <- .32 # rossington

##Load stan models from Stan models folder
load(here("Stan models","no_dist_seeds_pler_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_pler_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_pler_hi_lo.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_lo.Rdata"))
load(here("Stan models","no_dist_seeds_pler_lo_hi.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_hi.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_lo_hi.Rdata"))
load(here("Stan models","no_dist_seeds_pler_lo_int.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_int.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_lo_int.Rdata"))
load(here("Stan models","no_dist_seeds_pler_lo_lo.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_lo.Rdata"))

#Extract parameters from models and rename 
brho_hi_hi <- rstan::extract(no_dist_seeds_brho_hi_hi)
brho_hi_int <- rstan::extract(no_dist_seeds_brho_hi_int)
brho_hi_lo <- rstan::extract(no_dist_seeds_brho_hi_lo)
brho_lo_hi <- rstan::extract(no_dist_seeds_brho_lo_hi)
brho_lo_int <- rstan::extract(no_dist_seeds_brho_lo_int)
brho_lo_lo <- rstan::extract(no_dist_seeds_brho_lo_lo)
lapl_hi_hi <- rstan::extract(no_dist_seeds_lapl_hi_hi)
lapl_hi_int <- rstan::extract(no_dist_seeds_lapl_hi_int)
lapl_lo_hi <- rstan::extract(no_dist_seeds_lapl_lo_hi)
lapl_lo_int <- rstan::extract(no_dist_seeds_lapl_lo_int)
pler_hi_hi <- rstan::extract(no_dist_seeds_pler_hi_hi)
pler_hi_int <- rstan::extract(no_dist_seeds_pler_hi_int)
pler_hi_lo <- rstan::extract(no_dist_seeds_pler_hi_lo)
pler_lo_hi <- rstan::extract(no_dist_seeds_pler_lo_hi)
pler_lo_int <- rstan::extract(no_dist_seeds_pler_lo_int)
pler_lo_lo <- rstan::extract(no_dist_seeds_pler_lo_lo)

pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0*g)
  return(N)
}

################
######BRHO######
################

# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
lambdas <- brho_hi_hi$lambda[posts]
alphas_intra <- brho_hi_hi$alpha_brho[posts]
alphas_inter_pler <- brho_hi_hi$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_hi$alpha_lapl[posts]


## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = bs, g = bg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}

N_equil
N_brho_hi_hi <- N_equil[51]


##################
#######PLER#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
lambdas <- pler_hi_hi$lambda[posts]
alphas_intra <- pler_hi_hi$alpha_pler[posts]
alphas_inter_brho <- pler_hi_hi$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_hi$alpha_lapl[posts]


## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ps, g = pg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}
N_equil
N_pler_hi_hi <- N_equil[51]

##################
#######LAPL#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(lapl_hi_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
lambdas <- lapl_hi_hi$lambda[posts]
alphas_intra <- lapl_hi_hi$alpha_lapl[posts]
alphas_inter_brho <- lapl_hi_hi$alpha_brho[posts]
alphas_inter_pler <- lapl_hi_hi$alpha_pler[posts]

## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ls, g = lg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}
N_equil
N_lapl_hi_hi <- N_equil[51]

###########################################
####Equilibrium abundance HI H20 INT N#####
###########################################

################
######BRHO######
################

# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_int$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
lambdas <- brho_hi_int$lambda[posts]
alphas_intra <- brho_hi_int$alpha_brho[posts]
alphas_inter_pler <- brho_hi_int$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_int$alpha_lapl[posts]

## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = bs, g = bg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}

N_equil
N_brho_hi_int <- N_equil[51]


##################
#######PLER#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_int$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
lambdas <- pler_hi_int$lambda[posts]
alphas_intra <- pler_hi_int$alpha_pler[posts]
alphas_inter_brho <- pler_hi_int$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_int$alpha_lapl[posts]


## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ps, g = pg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}

N_equil
N_pler_hi_int <- N_equil[51]

##################
#######LAPL#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(lapl_hi_int$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
lambdas <- lapl_hi_int$lambda[posts]
alphas_intra <- lapl_hi_int$alpha_lapl[posts]
alphas_inter_brho <- lapl_hi_int$alpha_brho[posts]
alphas_inter_pler <- lapl_hi_int$alpha_pler[posts]

## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ls, g = lg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}
N_equil
N_lapl_hi_int <- N_equil[51]

###########################################
####Equilibrium abundance HI H20 LO N#####
###########################################

################
######BRHO######
################

# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
lambdas <- brho_hi_lo$lambda[posts]
alphas_intra <- brho_hi_lo$alpha_brho[posts]
alphas_inter_pler <- brho_hi_lo$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_lo$alpha_lapl[posts]

## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = bs, g = bg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}

# Calculate SE 
N_equil
N_brho_hi_lo <- N_equil[51]


##################
#######PLER#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
lambdas <- pler_hi_lo$lambda[posts]
alphas_intra <- pler_hi_lo$alpha_pler[posts]
alphas_inter_brho <- pler_hi_lo$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_lo$alpha_lapl[posts]


## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ps, g = pg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}

N_equil
N_pler_hi_lo <-  N_equil[51]

##################
#######LAPL#######
##################
N_lapl_hi_lo <- 0

##########################################
####Equilibrium abundance LO H20 HI N#####
##########################################

################
######BRHO######
################

# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
lambdas <- brho_lo_hi$lambda[posts]
alphas_intra <- brho_lo_hi$alpha_brho[posts]
alphas_inter_pler <- brho_lo_hi$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_hi$alpha_lapl[posts]


## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = bs, g = bg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}

# Calculate SE 
N_equil
N_brho_lo_hi <- N_equil[51]


##################
#######PLER#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
lambdas <- pler_lo_hi$lambda[posts]
alphas_intra <- pler_lo_hi$alpha_pler[posts]
alphas_inter_brho <- pler_lo_hi$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_hi$alpha_lapl[posts]


## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ps, g = pg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}
N_equil
N_pler_lo_hi <- N_equil[51]

##################
#######LAPL#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(lapl_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
lambdas <- lapl_lo_hi$lambda[posts]
alphas_intra <- lapl_lo_hi$alpha_lapl[posts]
alphas_inter_brho <- lapl_lo_hi$alpha_brho[posts]
alphas_inter_pler <- lapl_lo_hi$alpha_pler[posts]

## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ls, g = lg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}
N_equil
N_lapl_lo_hi <- N_equil[51]

###########################################
####Equilibrium abundance LO H20 INT N#####
###########################################

################
######BRHO######
################

# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_int$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
lambdas <- brho_lo_int$lambda[posts]
alphas_intra <- brho_lo_int$alpha_brho[posts]
alphas_inter_pler <- brho_lo_int$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_int$alpha_lapl[posts]

## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = bs, g = bg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}


N_equil
N_brho_lo_int <- N_equil[51]


##################
#######PLER#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_int$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
lambdas <- pler_lo_int$lambda[posts]
alphas_intra <- pler_lo_int$alpha_pler[posts]
alphas_inter_brho <- pler_lo_int$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_int$alpha_lapl[posts]


## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ps, g = pg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}

N_equil
N_pler_lo_int <- N_equil[51]

##################
#######LAPL#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(lapl_lo_int$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
lambdas <- lapl_lo_int$lambda[posts]
alphas_intra <- lapl_lo_int$alpha_lapl[posts]
alphas_inter_brho <- lapl_lo_int$alpha_brho[posts]
alphas_inter_pler <- lapl_lo_int$alpha_pler[posts]

## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ls, g = lg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}
N_equil
N_lapl_lo_int <- N_equil[51]

###########################################
####Equilibrium abundance LO H20 LO N#####
###########################################

################
######BRHO######
################

# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
lambdas <- brho_lo_lo$lambda[posts]
alphas_intra <- brho_lo_lo$alpha_brho[posts]
alphas_inter_pler <- brho_lo_lo$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_lo$alpha_lapl[posts]

## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 900

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = bs, g = bg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}


N_equil
N_brho_lo_lo <- N_equil[51]


##################
#######PLER#######
##################

# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
lambdas <- pler_lo_lo$lambda[posts]
alphas_intra <- pler_lo_lo$alpha_pler[posts]
alphas_inter_brho <- pler_lo_lo$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_lo$alpha_lapl[posts]


## For equilibrium abundance confidence intervals
# Set length of time for achieving equilibrium population
time <- 50

# Create output equilibrium vector
equil_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)) {
  
  # Create vector of abundances to reach equilibrium
  N_equil <- vector(length = time + 1)
  # Set starting population
  N_equil[1] <- 100
  
  # Loop once for each time step, ideally reaching equilibrium by final time step
  for (j in 1:time) {
    
    # Be sure to subset the ith values of posterior alphas and lambdas 
    N_equil[j + 1] <- pop.equilibrium(N0 = N_equil[j], s = ps, g = pg, 
                                      a_intra = alphas_intra[i], lambda = lambdas[i])
  }
  
  # Store final abundance in your output object
  equil_out[i] <- N_equil[j + 1]
  
}

N_equil
N_pler_lo_lo <- N_equil[51]

##################
#######LAPL#######
##################
N_lapl_lo_lo <- 0


#####put together into a dataframe
species_eq <- data.frame(N_brho_hi_hi,N_brho_hi_int,N_brho_hi_lo,N_brho_lo_hi,N_brho_lo_int,N_brho_lo_lo,
      N_lapl_hi_hi,N_lapl_hi_int,N_lapl_hi_lo,N_lapl_lo_hi,N_lapl_lo_int,N_lapl_lo_lo,
      N_pler_hi_hi,N_pler_hi_int,N_pler_hi_lo,N_pler_lo_hi,N_pler_lo_int,N_pler_lo_lo) %>%
  pivot_longer(1:18,names_to="trt",values_to = "equil_abundance") %>%
  separate(1,into=c("del","species","w_trt","n_trt"),sep="_") %>%
  select(-1) %>%
  mutate_at(4,as.integer)

species_eq$species[species_eq$species == "brho"] <- "Bromus"
species_eq$species[species_eq$species == "lapl"] <- "Layia"
species_eq$species[species_eq$species == "pler"] <- "Plantago"

species_eq$w_trt[species_eq$w_trt == "hi"] <- "Wet"
species_eq$w_trt[species_eq$w_trt == "lo"] <- "Dry"

species_eq$n_trt[species_eq$n_trt == "hi"] <- "hi.N"
species_eq$n_trt[species_eq$n_trt == "int"] <- "int.N"
species_eq$n_trt[species_eq$n_trt == "lo"] <- "lo.N"




