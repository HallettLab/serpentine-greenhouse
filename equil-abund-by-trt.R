##########################################
####Equilibrium abundance HI H20 Hi N#####
##########################################

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




