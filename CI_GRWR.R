#################
###GRWR Error####
#################

##Data
params <- read.csv(paste(datpath, "params3.csv", sep = ""))

##Functions##
#calculate SE
calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

# Determine equilibrium conditions for each species in isolation 
pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0)
  return(N)
}

# invader population growth rate one time step forward
pop.invade <- function (N0, resident, s, g, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_inter*resident)
  return(N)
}

# resident population growth rate one time step forward
pop.resident <- function (N0, resident, s, g, a_intra, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*resident + resident*(lambda*g)/(1+a_intra*resident+a_inter*N0)
  return(N)
}

##surv and germ fractions from literature
ps <- .75 # gulmon
pg <- .92 # gulmon
bs <- .013 # andrew
bg <- .98 # gulmon
ls <- .15 # ms thesis from Cal Poly SLO
lg <- .32 # ms thesis from Cal Poly SLO

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

# Calculate SE 
N_se <- calcSE(equil_out)
N_equil
N_brho_hi_hi <- 203.9877


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
N_pler_hi_hi <- 101.8302

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
N_lapl_hi_hi <- 119.5346

###########################
######CIs HI H2O HI N######
###########################
N_brho_hi_hi <- 203.9877
N_lapl_hi_hi <- 119.5346
N_pler_hi_hi <- 101.8302

##########################
###PLER invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_hi$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                             resident = N_brho_hi_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_brho_hi_hi_se <- ldgr_se

##########################
###PLER invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_hi_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_lapl_hi_hi_se <- ldgr_se

##########################
###LAPL invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(lapl_hi_hi$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_hi_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
lapl_brho_hi_hi_se <- ldgr_se

##########################
###LAPL invading PLER#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_hi_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
lapl_pler_hi_hi_se <- ldgr_se

##########################
###BRHO invading PLER#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_hi$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_hi_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_pler_hi_hi_se <- ldgr_se

##########################
###BRHO invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_hi_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_lapl_hi_hi_se <- ldgr_se

########################################################################################
########################################################################################

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

# Calculate SE 
N_se <- calcSE(equil_out)
N_brho_hi_int <- 61.63854


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
N_pler_hi_int <- 52.75433

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
N_lapl_hi_int <- 13.19490

###########################
######CIs HI H2O INT N######
###########################
N_brho_hi_int <- 61.63854
N_lapl_hi_int <- 13.19490
N_pler_hi_int <- 52.75433

##########################
###PLER invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_int$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_hi_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_brho_hi_int_se <- ldgr_se

##########################
###PLER invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_hi_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_lapl_hi_int_se <- ldgr_se

##########################
###LAPL invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(lapl_hi_int$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_hi_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
lapl_brho_hi_int_se <- ldgr_se

##########################
###LAPL invading PLER#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_hi_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
lapl_pler_hi_int_se <- ldgr_se

##########################
###BRHO invading PLER#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_hi$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_hi_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_pler_hi_int_se <- ldgr_se

##########################
###BRHO invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_hi_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_lapl_hi_int_se <- ldgr_se

########################################################################################
########################################################################################

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
N_se <- calcSE(equil_out)
N_equil
N_brho_hi_lo <- 15.57822


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
N_pler_hi_lo <- 60.52950

##################
#######LAPL#######
##################
N_lapl_hi_lo <- 0

###########################
######CIs HI H2O INT N######
###########################
N_brho_hi_lo <- 15.57822
N_lapl_hi_lo <- 0
N_pler_hi_lo <- 60.52950

##########################
###PLER invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_hi_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_brho_hi_lo_se <- ldgr_se

##########################
###PLER invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_hi_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_lapl_hi_lo_se <- ldgr_se


##########################
###BRHO invading PLER#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_hi_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_pler_hi_lo_se <- ldgr_se

##########################
###BRHO invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_hi_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_lapl_hi_lo_se <- ldgr_se

#####################################################################################
#####################################################################################

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
N_se <- calcSE(equil_out)
N_equil
N_brho_lo_hi <- 119.2125


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
N_pler_lo_hi <- 166.7332

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
N_lapl_lo_hi <- 122.6577

###########################
######CIs LO H2O HI N######
###########################
N_brho_lo_hi <- 119.2125
N_lapl_lo_hi <- 122.6577
N_pler_lo_hi <- 166.7332

##########################
###PLER invading BRHO#####
##########################

# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_lo_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_brho_lo_hi_se <- ldgr_se

##########################
###PLER invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_lo_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_lapl_lo_hi_se <- ldgr_se

##########################
###LAPL invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(lapl_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_lo_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
lapl_brho_lo_hi_se <- ldgr_se

##########################
###LAPL invading PLER#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_lo_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
lapl_pler_lo_hi_se <- ldgr_se

##########################
###BRHO invading PLER#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_hi$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_lo_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_pler_lo_hi_se <- ldgr_se

##########################
###BRHO invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_lo_hi, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_lapl_lo_hi_se <- ldgr_se

########################################################################################
########################################################################################

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

# Calculate SE 
N_se <- calcSE(equil_out)
N_equil
N_brho_lo_int <- 53.04479


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
N_pler_lo_int <- 57.82130

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
N_lapl_lo_int <- 4.434136

###########################
######CIs LO H2O INT N######
###########################
N_brho_lo_int <- 53.04479
N_lapl_lo_int <- 4.434136
N_pler_lo_int <- 57.82130

##########################
###PLER invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_int$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_lo_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_brho_lo_int_se <- ldgr_se

##########################
###PLER invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_lo_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_lapl_lo_int_se <- ldgr_se

##########################
###LAPL invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(lapl_lo_int$lambda)

# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)

# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_lo_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
lapl_brho_lo_int_se <- ldgr_se

##########################
###LAPL invading PLER#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ls, g = lg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_lo_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
lapl_pler_lo_int_se <- ldgr_se

##########################
###BRHO invading PLER#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_int$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_lo_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_pler_lo_int_se <- ldgr_se

##########################
###BRHO invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_lo_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_lapl_lo_int_se <- ldgr_se

########################################################################################
########################################################################################

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
N_se <- calcSE(equil_out)
N_equil
N_brho_hi_lo <- 15.57822


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
N_pler_hi_lo <- 60.52950

##################
#######LAPL#######
##################
N_lapl_hi_lo <- 0

###########################
######CIs HI H2O INT N######
###########################
N_brho_hi_lo <- 15.57822
N_lapl_hi_lo <- 0
N_pler_hi_lo <- 60.52950

##########################
###PLER invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(pler_hi_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_hi_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_brho_hi_lo_se <- ldgr_se

##########################
###PLER invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_hi_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_lapl_hi_lo_se <- ldgr_se


##########################
###BRHO invading PLER#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_hi_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_pler_hi_lo_se <- ldgr_se

##########################
###BRHO invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_hi_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_lapl_hi_lo_se <- ldgr_se

########################################################################################
########################################################################################

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

# Calculate SE 
N_se <- calcSE(equil_out)
N_equil
N_brho_lo_lo <- 0


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
N_pler_lo_lo <- 15.98185

##################
#######LAPL#######
##################
N_lapl_lo_lo <- 0

###########################
######CIs HI H2O INT N######
###########################
N_brho_lo_lo <- 0
N_lapl_lo_lo <- 0
N_pler_lo_lo <- 15.98185

##########################
###PLER invading BRHO#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(pler_lo_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_brho[i],
                            resident = N_brho_lo_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_brho_lo_lo_se <- ldgr_se

##########################
###PLER invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = ps, g = pg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_lo_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
pler_lapl_lo_lo_se <- ldgr_se


##########################
###BRHO invading PLER#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(brho_lo_lo$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_lo_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_pler_lo_lo_se <- ldgr_se

##########################
###BRHO invading LAPL#####
##########################
## For LDGR confidence intervals
# Create output object for invader LDGRs
ldgr_out <- vector(length = length(posts))

# Loop for once set of randomly sampled parameters from your poteriors
for (i in 1:length(posts)){
  
  # Calculate an LDGR for each set of sampled lambdas and alphas
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_lapl[i],
                            resident = N_lapl_lo_lo, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- calcSE(ldgr_out)
brho_lapl_lo_lo_se <- ldgr_se

grwr_se <- data.frame(brho_lapl_hi_hi_se,brho_lapl_hi_int_se,brho_lapl_hi_lo_se,
                      brho_lapl_lo_hi_se,brho_lapl_lo_int_se,brho_lapl_lo_lo_se,
                      pler_brho_hi_hi_se,pler_brho_hi_int_se,pler_brho_hi_lo_se,
                      pler_brho_lo_hi_se,pler_brho_lo_int_se,pler_brho_lo_lo_se,
                      lapl_brho_hi_hi_se,lapl_brho_hi_int_se,lapl_brho_lo_hi_se,
                      lapl_brho_lo_int_se,brho_pler_hi_hi_se,brho_pler_hi_int_se,
                      brho_pler_hi_lo_se,brho_pler_lo_hi_se,brho_pler_lo_int_se,
                      brho_pler_lo_lo_se,pler_lapl_hi_hi_se,pler_lapl_hi_int_se,
                      pler_lapl_hi_lo_se,pler_lapl_lo_hi_se,pler_lapl_lo_int_se,
                      pler_lapl_lo_lo_se,lapl_pler_hi_hi_se,lapl_pler_hi_int_se,
                      lapl_pler_lo_hi_se,lapl_pler_lo_int_se)

grwr_se<-pivot_longer(grwr_se,cols=1:32,names_to = "treatment", values_to = "se")
