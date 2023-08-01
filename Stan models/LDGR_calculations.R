# Equations and literature fractions ----
# example script for calculating LDGR with the posteriors

library(tidyverse)
library(rstan)
library(here)

# Determine equilibrium conditions for each species in isolation 
# LGS: Note I added *g in the denominator to switch from seeds NO to stems (NO*g)
pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0*g)
  return(N)
}

# invader population growth rate one time step forward
# LGS: Note I added *g_resident in the denominator to switch from seeds NO 
# to stems (NO*g). This also then means I addedd g_resident as an argument
# to the function
pop.invade <- function (N0, resident, s, g, a_inter, lambda, g_resident) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_inter*resident*g_resident)
  return(N)
}

##surv and germ fractions from literature
ps <- .75 # gulmon
pg <- .92 # gulmon
bs <- .013 # andrew
bg <- .98 # gulmon
ls <- .15 # ms thesis from Cal Poly SLO
lg <- .32 # ms thesis from Cal Poly SLO

# High water high N LDGR ----

# extract posteriors
load(here("Stan models","no_dist_seeds_pler_hi_hi.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_hi.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_hi_hi.Rdata"))
pler_hi_hi <- rstan::extract(no_dist_seeds_pler_hi_hi)
brho_hi_hi <- rstan::extract(no_dist_seeds_brho_hi_hi)
lapl_hi_hi <- rstan::extract(no_dist_seeds_lapl_hi_hi)


# BRHO to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_hi$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_brho_hi_hi <-equil_out

# PLER to equilibrium 

# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_hi$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_pler_hi_hi <-equil_out

# LAPL to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(lapl_hi_hi$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_lapl_hi_hi <-equil_out

#PLER invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_hi_hi$lambda[posts]
alphas_intra <- pler_hi_hi$alpha_pler[posts]
alphas_inter_brho <- pler_hi_hi$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_hi_hi[i], N0 = 1, g_resident=bg ))
  
}

pler_brho_hi_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_brho_hi_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of pler invading brho at hi_hi conditions.

#BRHO invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_hi_hi$lambda[posts]
alphas_intra <- brho_hi_hi$alpha_brho[posts]
alphas_inter_pler <- brho_hi_hi$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_hi_hi[i], N0 = 1, g_resident=pg ))
  
}

brho_pler_hi_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
brho_pler_hi_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading pler at hi_hi conditions.

#LAPL invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- lapl_hi_hi$lambda[posts]
alphas_intra <- lapl_hi_hi$alpha_lapl[posts]
alphas_inter_brho <- lapl_hi_hi$alpha_brho[posts]
alphas_inter_pler <- lapl_hi_hi$alpha_pler[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_hi_hi[i], N0 = 1, g_resident=bg ))
  
}

lapl_brho_hi_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_brho_hi_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

#BRHO invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_hi_hi$lambda[posts]
alphas_intra <- brho_hi_hi$alpha_brho[posts]
alphas_inter_pler <- brho_hi_hi$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_hi_hi[i], N0 = 1, g_resident=lg ))
  
}

brho_lapl_hi_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
brho_lapl_hi_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading lapl at hi_hi conditions.

#LAPL invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- lapl_hi_hi$lambda[posts]
alphas_intra <- lapl_hi_hi$alpha_lapl[posts]
alphas_inter_brho <- lapl_hi_hi$alpha_brho[posts]
alphas_inter_pler <- lapl_hi_hi$alpha_pler[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_hi_hi[i], N0 = 1, g_resident=pg ))
  
}

lapl_pler_hi_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_pler_hi_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

#PLER invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_hi_hi$lambda[posts]
alphas_intra <- pler_hi_hi$alpha_pler[posts]
alphas_inter_brho <- pler_hi_hi$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_hi_hi[i], N0 = 1, g_resident=lg ))
  
}

pler_lapl_hi_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_lapl_hi_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.


# High water intermediate N LDGR ----

# extract posteriors
load(here("Stan models","no_dist_seeds_pler_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_hi_int.Rdata"))
pler_hi_int <- rstan::extract(no_dist_seeds_pler_hi_int)
brho_hi_int <- rstan::extract(no_dist_seeds_brho_hi_int)
lapl_hi_int <- rstan::extract(no_dist_seeds_lapl_hi_int)


# BRHO to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_int$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_brho_hi_int <-equil_out

# PLER to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_int$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_pler_hi_int <-equil_out

#LAPL to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(lapl_hi_int$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_lapl_hi_int <-equil_out

# PLER invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_hi_int$lambda[posts]
alphas_intra <- pler_hi_int$alpha_pler[posts]
alphas_inter_brho <- pler_hi_int$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_hi_int[i], N0 = 1, g_resident=bg ))
  
}

pler_brho_hi_int_sd <- sd(ldgr_out,na.rm=TRUE) # LGS: since we take the log above, we are all set here
pler_brho_hi_int_mean <- mean(ldgr_out,na.rm=TRUE) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of pler invading brho at hi_hi conditions.

# BRHO invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_hi_int$lambda[posts]
alphas_intra <- brho_hi_int$alpha_brho[posts]
alphas_inter_pler <- brho_hi_int$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_int$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_hi_int[i], N0 = 1, g_resident=pg ))
  
}

brho_pler_hi_int_sd <- sd(ldgr_out,na.rm=TRUE) # LGS: since we take the log above, we are all set here
brho_pler_hi_int_mean <- mean(ldgr_out,na.rm=TRUE) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading pler at hi_hi conditions.

# LAPL invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- lapl_hi_int$lambda[posts]
alphas_intra <- lapl_hi_int$alpha_lapl[posts]
alphas_inter_brho <- lapl_hi_int$alpha_brho[posts]
alphas_inter_pler <- lapl_hi_int$alpha_pler[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_hi_int[i], N0 = 1, g_resident=bg ))
  
}

lapl_brho_hi_int_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_brho_hi_int_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# BRHO invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_hi_int$lambda[posts]
alphas_intra <- brho_hi_int$alpha_brho[posts]
alphas_inter_pler <- brho_hi_int$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_int$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_hi_int[i], N0 = 1, g_resident=lg ))
  
}

brho_lapl_hi_int_sd <- sd(ldgr_out,na.rm=TRUE) # LGS: since we take the log above, we are all set here
brho_lapl_hi_int_mean <- mean(ldgr_out,na.rm=TRUE) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading lapl at hi_hi conditions.

# LAPL invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- lapl_hi_int$lambda[posts]
alphas_intra <- lapl_hi_int$alpha_lapl[posts]
alphas_inter_brho <- lapl_hi_int$alpha_brho[posts]
alphas_inter_pler <- lapl_hi_int$alpha_pler[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_hi_int[i], N0 = 1, g_resident=pg ))
  
}

lapl_pler_hi_int_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_pler_hi_int_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# PLER invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_hi_int$lambda[posts]
alphas_intra <- pler_hi_int$alpha_pler[posts]
alphas_inter_brho <- pler_hi_int$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_int$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_hi_int[i], N0 = 1, g_resident=lg ))
  
}

pler_lapl_hi_int_sd <- sd(ldgr_out,na.rm=TRUE) # LGS: since we take the log above, we are all set here
pler_lapl_hi_int_mean <- mean(ldgr_out,na.rm=TRUE) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# High water low N LDGR ----

# extract posteriors
load(here("Stan models","no_dist_seeds_pler_hi_lo.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_lo.Rdata"))
pler_hi_lo <- rstan::extract(no_dist_seeds_pler_hi_lo)
brho_hi_lo <- rstan::extract(no_dist_seeds_brho_hi_lo)

# BRHO to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_lo$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_brho_hi_lo <-equil_out

# PLER to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_lo$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_pler_hi_lo <-equil_out

# LAPL to equilibrium

N_lapl_hi_lo <- 0

# PLER invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_hi_lo$lambda[posts]
alphas_intra <- pler_hi_lo$alpha_pler[posts]
alphas_inter_brho <- pler_hi_lo$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_lo$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_hi_lo[i], N0 = 1, g_resident=bg ))
  
}

pler_brho_hi_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_brho_hi_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of pler invading brho at hi_hi conditions.

# BRHO invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_hi_lo$lambda[posts]
alphas_intra <- brho_hi_lo$alpha_brho[posts]
alphas_inter_pler <- brho_hi_lo$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_lo$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_hi_lo[i], N0 = 1, g_resident=pg ))
  
}

brho_pler_hi_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
brho_pler_hi_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading pler at hi_hi conditions.

# LAPL invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- 0
alphas_intra <- 0
alphas_inter_brho <- 0
alphas_inter_pler <- 0

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas, a_inter = alphas_inter_brho,
                                resident = N_brho_hi_lo[i], N0 = 1, g_resident=bg ))
  
}

lapl_brho_hi_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_brho_hi_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# BRHO invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_hi_lo$lambda[posts]
alphas_intra <- brho_hi_lo$alpha_brho[posts]
alphas_inter_pler <- brho_hi_lo$alpha_pler[posts]
alphas_inter_lapl <- brho_hi_lo$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_hi_lo, N0 = 1, g_resident=lg ))
  
}

brho_lapl_hi_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
brho_lapl_hi_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading lapl at hi_hi conditions.

# LAPL invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- 0
alphas_intra <- 0
alphas_inter_brho <- 0
alphas_inter_pler <- 0

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas, a_inter = alphas_inter_pler,
                                resident = N_pler_hi_lo[i], N0 = 1, g_resident=pg ))
  
}

lapl_pler_hi_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_pler_hi_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# PLER invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_hi_lo$lambda[posts]
alphas_intra <- pler_hi_lo$alpha_pler[posts]
alphas_inter_brho <- pler_hi_lo$alpha_brho[posts]
alphas_inter_lapl <- pler_hi_lo$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_hi_lo, N0 = 1, g_resident=lg ))
  
}

pler_lapl_hi_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_lapl_hi_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# Low water high N LDGR ----

# extract posteriors
load(here("Stan models","no_dist_seeds_pler_lo_hi.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_hi.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_lo_hi.Rdata"))
pler_lo_hi <- rstan::extract(no_dist_seeds_pler_lo_hi)
brho_lo_hi <- rstan::extract(no_dist_seeds_brho_lo_hi)
lapl_lo_hi <- rstan::extract(no_dist_seeds_lapl_lo_hi)

# BRHO to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_brho_lo_hi <-equil_out

# PLER to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_pler_lo_hi <-equil_out

# LAPL to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(lapl_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_lapl_lo_hi <-equil_out

# PLER invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_lo_hi$lambda[posts]
alphas_intra <- pler_lo_hi$alpha_pler[posts]
alphas_inter_brho <- pler_lo_hi$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_lo_hi[i], N0 = 1, g_resident=bg ))
  
}

pler_brho_lo_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_brho_lo_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of pler invading brho at hi_hi conditions.

# BRHO invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_lo_hi$lambda[posts]
alphas_intra <- brho_lo_hi$alpha_brho[posts]
alphas_inter_pler <- brho_lo_hi$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_lo_hi[i], N0 = 1, g_resident=pg ))
  
}

brho_pler_lo_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
brho_pler_lo_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading pler at hi_hi conditions.

# LAPL invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- lapl_lo_hi$lambda[posts]
alphas_intra <- lapl_lo_hi$alpha_lapl[posts]
alphas_inter_brho <- lapl_lo_hi$alpha_brho[posts]
alphas_inter_pler <- lapl_lo_hi$alpha_pler[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_lo_hi[i], N0 = 1, g_resident=bg ))
  
}

lapl_brho_lo_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_brho_lo_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# BRHO invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_lo_hi$lambda[posts]
alphas_intra <- brho_lo_hi$alpha_brho[posts]
alphas_inter_pler <- brho_lo_hi$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_lo_hi[i], N0 = 1, g_resident=lg ))
  
}

brho_lapl_lo_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
brho_lapl_lo_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading lapl at hi_hi conditions.

# LAPL invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- lapl_lo_hi$lambda[posts]
alphas_intra <- lapl_lo_hi$alpha_lapl[posts]
alphas_inter_brho <- lapl_lo_hi$alpha_brho[posts]
alphas_inter_pler <- lapl_lo_hi$alpha_pler[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_lo_hi[i], N0 = 1, g_resident=pg ))
  
}

lapl_pler_lo_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_pler_lo_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# PLER invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_lo_hi$lambda[posts]
alphas_intra <- pler_lo_hi$alpha_pler[posts]
alphas_inter_brho <- pler_lo_hi$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_hi$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_lo_hi[i], N0 = 1, g_resident=lg ))
  
}

pler_lapl_lo_hi_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_lapl_lo_hi_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# Low water intermediate N LDGR ----
# extract posteriors
load(here("Stan models","no_dist_seeds_pler_lo_int.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_int.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_lo_int.Rdata"))
pler_lo_int <- rstan::extract(no_dist_seeds_pler_lo_int)
brho_lo_int <- rstan::extract(no_dist_seeds_brho_lo_int)
lapl_lo_int <- rstan::extract(no_dist_seeds_lapl_lo_int)


#BRHO to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_int$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_brho_lo_int <-equil_out

# PLER to euilibrium

# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_int$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_pler_lo_int <-equil_out

# LAPL to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(lapl_lo_int$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_lapl_lo_int <-equil_out

# PLER invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_lo_int$lambda[posts]
alphas_intra <- pler_lo_int$alpha_pler[posts]
alphas_inter_brho <- pler_lo_int$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_int$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_lo_int[i], N0 = 1, g_resident=bg ))
  
}

pler_brho_lo_int_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_brho_lo_int_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of pler invading brho at hi_hi conditions.

# BRHO invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_lo_int$lambda[posts]
alphas_intra <- brho_lo_int$alpha_brho[posts]
alphas_inter_pler <- brho_lo_int$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_int$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_lo_int[i], N0 = 1, g_resident=pg ))
  
}

brho_pler_lo_int_sd <- sd(ldgr_out,na.rm=TRUE) # LGS: since we take the log above, we are all set here
brho_pler_lo_int_mean <- mean(ldgr_out,na.rm=TRUE) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading pler at hi_hi conditions.

# LAPL invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- lapl_lo_int$lambda[posts]
alphas_intra <- lapl_lo_int$alpha_lapl[posts]
alphas_inter_brho <- lapl_lo_int$alpha_brho[posts]
alphas_inter_pler <- lapl_lo_int$alpha_pler[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_lo_int[i], N0 = 1, g_resident=bg ))
  
}

lapl_brho_lo_int_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_brho_lo_int_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# BRHO invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_lo_int$lambda[posts]
alphas_intra <- brho_lo_int$alpha_brho[posts]
alphas_inter_pler <- brho_lo_int$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_int$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_lo_int[i], N0 = 1, g_resident=lg ))
  
}

brho_lapl_lo_int_sd <- sd(ldgr_out,na.rm=TRUE) # LGS: since we take the log above, we are all set here
brho_lapl_lo_int_mean <- mean(ldgr_out,na.rm=TRUE) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading lapl at hi_hi conditions.

# LAPL invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- lapl_lo_int$lambda[posts]
alphas_intra <- lapl_lo_int$alpha_lapl[posts]
alphas_inter_brho <- lapl_lo_int$alpha_brho[posts]
alphas_inter_pler <- lapl_lo_int$alpha_pler[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_lo_int[i], N0 = 1, g_resident=pg ))
  
}

lapl_pler_lo_int_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_pler_lo_int_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# PLER invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_lo_int$lambda[posts]
alphas_intra <- pler_lo_int$alpha_pler[posts]
alphas_inter_brho <- pler_lo_int$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_int$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_lo_int[i], N0 = 1, g_resident=lg ))
  
}

pler_lapl_lo_int_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_lapl_lo_int_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# Low water low N LDGR ----

# extract posteriors
load(here("Stan models","no_dist_seeds_pler_lo_lo.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_lo.Rdata"))
pler_lo_lo <- rstan::extract(no_dist_seeds_pler_lo_lo)
brho_lo_lo <- rstan::extract(no_dist_seeds_brho_lo_lo)


# BRHO to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_lo$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

# Subset your posteriors using your random positions
lambdas <- brho_lo_lo$lambda[posts]
alphas_intra <- brho_lo_lo$alpha_brho[posts]
alphas_inter_pler <- brho_lo_lo$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_lo$alpha_lapl[posts]


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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_brho_lo_lo <-equil_out

# PLER to equilibrium

# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_lo$lambda)

# Grab 400 random positions from your posterior distribution
# LGS let's increase this to 2000
# Also changed to replace = FALSE
posts <- sample(posterior_length, 2000, replace=FALSE)

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

# LGS: Let's save our resident equilibrium results, so we can use these
# for the invasion analyses below
N_pler_lo_lo <-equil_out

# LAPL to equilibrium
N_lapl_lo_lo <- 0

# PLER invading BRHO


# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- pler_lo_lo$lambda[posts]
alphas_intra <- pler_lo_lo$alpha_pler[posts]
alphas_inter_brho <- pler_lo_lo$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_lo$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                                resident = N_brho_lo_lo[i], N0 = 1, g_resident=bg ))
  
}

pler_brho_lo_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_brho_lo_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of pler invading brho at hi_hi conditions.

# BRHO invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_lo_lo$lambda[posts]
alphas_intra <- brho_lo_lo$alpha_brho[posts]
alphas_inter_pler <- brho_lo_lo$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_lo$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                                resident = N_pler_lo_lo[i], N0 = 1, g_resident=pg ))
  
}

brho_pler_lo_lo_sd <- sd(ldgr_out,na.rm=TRUE) # LGS: since we take the log above, we are all set here
brho_pler_lo_lo_mean <- mean(ldgr_out,na.rm=TRUE) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading pler at hi_hi conditions.

# LAPL invading BRHO

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- 0
alphas_intra <- 0
alphas_inter_brho <- 0
alphas_inter_pler <- 0

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas, a_inter = alphas_inter_brho,
                                resident = N_brho_lo_lo[i], N0 = 1, g_resident=bg ))
  
}

lapl_brho_lo_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_brho_lo_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# BRHO invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- brho_lo_lo$lambda[posts]
alphas_intra <- brho_lo_lo$alpha_brho[posts]
alphas_inter_pler <- brho_lo_lo$alpha_pler[posts]
alphas_inter_lapl <- brho_lo_lo$alpha_lapl[posts]

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_lo_lo, N0 = 1, g_resident=lg ))
  
}

brho_lapl_lo_lo_sd <- sd(ldgr_out,na.rm=TRUE) # LGS: since we take the log above, we are all set here
brho_lapl_lo_lo_mean <- mean(ldgr_out,na.rm=TRUE) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of brho invading lapl at hi_hi conditions.

# LAPL invading PLER

# Note example where both are at hi_hi
# We'll keep the same posts as above
lambdas <- 0
alphas_intra <- 0
alphas_inter_brho <- 0
alphas_inter_pler <- 0

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ls, g = lg, lambda=lambdas, a_inter = alphas_inter_pler,
                                resident = N_pler_lo_lo[i], N0 = 1, g_resident=pg ))
  
}

lapl_pler_lo_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
lapl_pler_lo_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.

# PLER invading LAPL

# Note example where both are at hi_hi
# We'll keep the same posts as above

## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))
lambdas <- pler_lo_lo$lambda[posts]
alphas_intra <- pler_lo_lo$alpha_pler[posts]
alphas_inter_brho <- pler_lo_lo$alpha_brho[posts]
alphas_inter_lapl <- pler_lo_lo$alpha_lapl[posts]

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # LGS: I just did the log calculation here
  # technically, it is log(Nt+1/Nt), but since Nt=1, we can just do log(Nt+1), as
  # we do here
  #LGS: note I added N_brho_hi_hi[i] with the i index to loop through our posteriors
  # of brho's equilibrium abundance.
  ldgr_out[i] <- log(pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                                resident = N_lapl_lo_lo, N0 = 1, g_resident=lg ))
  
}

pler_lapl_lo_lo_sd <- sd(ldgr_out) # LGS: since we take the log above, we are all set here
pler_lapl_lo_lo_mean <- mean(ldgr_out) # LGS: same comment as above

# LGS, this now gives us a posterior distribution (ldgr_out) of all expected
# ldgr of lapl invading brho at hi_hi conditions.


##############################################
##make dataframe for mean and sd GRWR values##
##############################################

mean_LDGR <- as.data.frame(cbind(brho_lapl_hi_hi_mean, brho_lapl_hi_int_mean, brho_lapl_hi_lo_mean,
                           brho_lapl_lo_hi_mean, brho_lapl_lo_int_mean, brho_lapl_lo_lo_mean,
                           lapl_brho_hi_hi_mean, lapl_brho_hi_int_mean, lapl_brho_hi_lo_mean,
                           lapl_brho_lo_hi_mean, lapl_brho_lo_int_mean, lapl_brho_lo_lo_mean,
                           brho_pler_hi_hi_mean, brho_pler_hi_int_mean, brho_pler_hi_lo_mean,
                           brho_pler_lo_hi_mean, brho_pler_lo_int_mean, brho_pler_lo_lo_mean,
                           pler_brho_hi_hi_mean, pler_brho_hi_int_mean, pler_brho_hi_lo_mean,
                           pler_brho_lo_hi_mean, pler_brho_lo_int_mean, pler_brho_lo_lo_mean,
                           lapl_pler_hi_hi_mean, lapl_pler_hi_int_mean, lapl_pler_hi_lo_mean,
                           lapl_pler_lo_hi_mean, lapl_pler_lo_int_mean, lapl_pler_lo_lo_mean,
                           pler_lapl_hi_hi_mean, pler_lapl_hi_int_mean, pler_lapl_hi_lo_mean,
                           pler_lapl_lo_hi_mean, pler_lapl_lo_int_mean, pler_lapl_lo_lo_mean)) %>%
  pivot_longer(1:36,names_to = "trt",values_to = "mean_LDGR") %>%
  separate(trt,into = c("invader","resident","water_trt","n_trt","del"),sep = "_")%>%
  unite("invader_resident",1:2) %>%
  select(-4)

sd_LDGR <- as.data.frame(cbind(brho_lapl_hi_hi_sd, brho_lapl_hi_int_sd, brho_lapl_hi_lo_sd,
                                 brho_lapl_lo_hi_sd, brho_lapl_lo_int_sd, brho_lapl_lo_lo_sd,
                                 lapl_brho_hi_hi_sd, lapl_brho_hi_int_sd, lapl_brho_hi_lo_sd,
                                 lapl_brho_lo_hi_sd, lapl_brho_lo_int_sd, lapl_brho_lo_lo_sd,
                                 brho_pler_hi_hi_sd, brho_pler_hi_int_sd, brho_pler_hi_lo_sd,
                                 brho_pler_lo_hi_sd, brho_pler_lo_int_sd, brho_pler_lo_lo_sd,
                                 pler_brho_hi_hi_sd, pler_brho_hi_int_sd, pler_brho_hi_lo_sd,
                                 pler_brho_lo_hi_sd, pler_brho_lo_int_sd, pler_brho_lo_lo_sd,
                                 lapl_pler_hi_hi_sd, lapl_pler_hi_int_sd, lapl_pler_hi_lo_sd,
                                 lapl_pler_lo_hi_sd, lapl_pler_lo_int_sd, lapl_pler_lo_lo_sd,
                                 pler_lapl_hi_hi_sd, pler_lapl_hi_int_sd, pler_lapl_hi_lo_sd,
                                 pler_lapl_lo_hi_sd, pler_lapl_lo_int_sd, pler_lapl_lo_lo_sd)) %>%
  pivot_longer(1:36,names_to = "trt",values_to = "sd_LDGR") %>%
  separate(trt,into = c("invader","resident","water_trt","n_trt","del"),sep = "_")%>%
  unite("invader_resident",1:2) %>%
  select(-4)

mean_sd_LDGR <- left_join(mean_LDGR,sd_LDGR)

###################################################
###LDGR/GRWR data manipulation & visualization#####
###################################################
LDGR <- mean_sd_LDGR %>%
  mutate(invader_resident2=invader_resident) %>%
  separate(1,into="invader",sep="_") %>%
  rename(invader_resident = invader_resident2) %>%
  unite("invader_water",1:2,sep = "_",remove = FALSE)

LDGR$invader_resident <- factor(LDGR$invader_resident, 
                                        levels = c("pler_brho","pler_lapl","lapl_brho",
                                                   "brho_pler","lapl_pler","brho_lapl"))
LDGR$n_trt <- factor(mean_sd_LDGR$n_trt, levels = c("lo","int","hi"))
LDGR$water_trt <- factor(mean_sd_LDGR$water_trt, levels = c("lo","hi"))

theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 16),
              axis.text = element_text(size = 16))

grwr <- ggplot(LDGR, aes(x=n_trt, y=mean_LDGR,ymin=mean_LDGR-sd_LDGR,
                                 ymax=mean_LDGR+sd_LDGR, fill = invader_water)) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  facet_wrap(~invader_resident) +
  theme(strip.text.x = element_blank(),legend.position= "none",
        axis.title = element_blank(),axis.text.x = element_blank())+
  ylab("Growth rate when rare") + xlab("N treatments") +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  scale_fill_manual(name="Water treatments", 
                    values=c("#D55E00","#f79f59","#0072B2","#49a0d1","#009E73","#5fc9ac")) +
  geom_hline(yintercept = 0)+
  geom_errorbar(position = position_dodge(0.9),width=0.1)

pdf("grwr.pdf", width = 7, height = 5)
grwr
dev.off()
