#################
###GRWR Error####
#################

##Data
params <- read.csv(paste(datpath, "params2.csv", sep = ""))

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

## extract from stan models
pler_lo_lo <- rstan::extract(no_dist_seeds_pler_lo_lo)
brho_lo_lo <- rstan::extract(no_dist_seeds_brho_lo_lo)
pler_lo_int <- rstan::extract(no_dist_seeds_pler_lo_int)
brho_lo_int <- rstan::extract(no_dist_seeds_brho_lo_int)
lapl_lo_int <- rstan::extract(no_dist_seeds_lapl_lo_int)
pler_lo_hi <- rstan::extract(no_dist_seeds_pler_lo_hi)
brho_lo_hi <- rstan::extract(no_dist_seeds_brho_lo_hi)
lapl_lo_hi <- rstan::extract(no_dist_seeds_lapl_lo_hi)
pler_hi_lo <- rstan::extract(no_dist_seeds_pler_hi_lo)
brho_hi_lo <- rstan::extract(no_dist_seeds_brho_hi_lo)
pler_hi_int <- rstan::extract(no_dist_seeds_pler_hi_int)
brho_hi_int <- rstan::extract(no_dist_seeds_brho_hi_int)
lapl_hi_int <- rstan::extract(no_dist_seeds_lapl_hi_int)
pler_hi_hi <- rstan::extract(no_dist_seeds_pler_hi_hi)
brho_hi_hi <- rstan::extract(no_dist_seeds_brho_hi_hi)
lapl_hi_hi <- rstan::extract(no_dist_seeds_lapl_hi_hi)

#two-step fit models
pler_lo_lo <- rstan::extract(pler_lo_lo)
brho_lo_lo <- rstan::extract(brho_lo_lo)
pler_lo_int <- rstan::extract(pler_lo_int)
brho_lo_int <- rstan::extract(brho_lo_int)
lapl_lo_int <- rstan::extract(lapl_lo_int)
pler_lo_hi <- rstan::extract(pler_lo_hi)
brho_lo_hi <- rstan::extract(brho_lo_hi)
lapl_lo_hi <- rstan::extract(lapl_lo_hi)
pler_hi_lo <- rstan::extract(pler_hi_lo)
brho_hi_lo <- rstan::extract(brho_hi_lo)
pler_hi_int <- rstan::extract(pler_hi_int)
brho_hi_int <- rstan::extract(brho_hi_int)
lapl_hi_int <- rstan::extract(lapl_hi_int)
pler_hi_hi <- rstan::extract(pler_hi_hi)
brho_hi_hi <- rstan::extract(brho_hi_hi)
lapl_hi_hi <- rstan::extract(lapl_hi_hi)

#pairs plots og stan model fits
#pler
pler_hi_hi<-data.frame(pler_hi_hi)
pler_hi_int<-data.frame(pler_hi_int)
pler_hi_lo<-data.frame(pler_hi_lo)
pler_lo_hi<-data.frame(pler_lo_hi)
pler_lo_int<-data.frame(pler_lo_int)
pler_lo_lo<-data.frame(pler_lo_lo)
plot_grid(
ggmatrix_gtable(ggpairs(pler_hi_hi[6:9],title="pler hi hi")),
  ggmatrix_gtable(ggpairs(pler_hi_int[6:9],title="pler hi int")),
  ggmatrix_gtable(ggpairs(pler_hi_lo[6:9],title="pler hi lo"),
  ggmatrix_gtable(ggpairs(pler_lo_hi[6:9],title="pler lo hi")),
  ggmatrix_gtable(ggpairs(pler_lo_int[6:9],title="pler lo int")),
  ggmatrix_gtable(ggpairs(pler_lo_lo[6:9],title= "pler lo lo")),
  nrow = 2,
  ncol = 3)

#brho
brho_hi_hi<-data.frame(brho_hi_hi)
brho_hi_int<-data.frame(brho_hi_int)
brho_hi_lo<-data.frame(brho_hi_lo)
brho_lo_hi<-data.frame(brho_lo_hi)
brho_lo_int<-data.frame(brho_lo_int)
brho_lo_lo<-data.frame(brho_lo_lo)
plot_grid(
  ggmatrix_gtable(ggpairs(brho_hi_hi[6:9],title="brho hi hi")),
  ggmatrix_gtable(ggpairs(brho_hi_int[6:9],title="brho hi int")),
  ggmatrix_gtable(ggpairs(brho_hi_lo[6:9],title="brho hi lo")),
  ggmatrix_gtable(ggpairs(brho_lo_hi[6:9],title="brho lo hi")),
  ggmatrix_gtable(ggpairs(brho_lo_int[6:9],title="brho lo int")),
  ggmatrix_gtable(ggpairs(brho_lo_lo[6:9],title="brho lo lo")),
  nrow = 2,
  ncol = 3)

#lapl
lapl_hi_hi<-data.frame(lapl_hi_hi)
lapl_hi_int<-data.frame(lapl_hi_int)
lapl_lo_hi<-data.frame(lapl_lo_hi)
lapl_lo_int<-data.frame(lapl_lo_int)
plot_grid(
  ggmatrix_gtable(ggpairs(lapl_hi_hi[6:9],title="lapl hi hi")),
  ggmatrix_gtable(ggpairs(lapl_hi_int[6:9],title="lapl hi int")),
  ggmatrix_gtable(ggpairs(lapl_lo_hi[6:9],title="lapl lo hi")),
  ggmatrix_gtable(ggpairs(brho_lo_int[6:9],title="lapl lo int")),
  nrow = 2,
  ncol = 2)


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
N_equil
N_brho_hi_hi <- 207.2087


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
N_pler_hi_hi <- 189.1819

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
N_lapl_hi_hi <- 140.6358


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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
pler_brho_hi_hi_se <- ldgr_se
pler_brho_hi_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
pler_lapl_hi_hi_se <- ldgr_se
pler_lapl_hi_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
lapl_brho_hi_hi_se <- ldgr_se
lapl_brho_hi_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
lapl_pler_hi_hi_se <- ldgr_se
lapl_pler_hi_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
brho_pler_hi_hi_se <- ldgr_se
brho_pler_hi_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
brho_lapl_hi_hi_se <- ldgr_se
brho_lapl_hi_hi_mean <- ldgr_mean

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
N_equil
N_brho_hi_int <- 57.11660


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
N_pler_hi_int <- 48.26259

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
N_lapl_hi_int <- 6.627693


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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
pler_brho_hi_int_se <- ldgr_se
pler_brho_hi_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
pler_lapl_hi_int_se <- ldgr_se
pler_lapl_hi_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
lapl_brho_hi_int_se <- ldgr_se
lapl_brho_hi_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
lapl_pler_hi_int_se <- ldgr_se
lapl_pler_hi_int_mean <- ldgr_mean

##########################
###BRHO invading PLER#####
##########################
# Grab the length of your posterior distribution
posterior_length <- length(brho_hi_int$lambda)
# Grab 400 random positions from your posterior distribution
posts <- sample(posterior_length, 400, replace=TRUE)
# Subset your posteriors using your random positions
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
  # Here your lambdas will always be the same, but the alpha_inter will change depending on the resident species
  ldgr_out[i] <- pop.invade(s = bs, g = bg, lambda=lambdas[i], a_inter = alphas_inter_pler[i],
                            resident = N_pler_hi_int, N0 = 1)
  
}

# Calculate SE 
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
brho_pler_hi_int_se <- ldgr_se
brho_pler_hi_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
brho_lapl_hi_int_se <- ldgr_se
brho_lapl_hi_int_mean <- ldgr_mean

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
N_equil
N_brho_hi_lo <- 17.93910


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
N_pler_hi_lo <-  72.57050

##################
#######LAPL#######
##################
N_lapl_hi_lo <- 0


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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
pler_brho_hi_lo_se <- ldgr_se
pler_brho_hi_lo_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
pler_lapl_hi_lo_se <- ldgr_se
pler_lapl_hi_lo_mean <- ldgr_mean


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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean(ldgr_out)
brho_pler_hi_lo_se <- ldgr_se
brho_pler_hi_lo_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
brho_lapl_hi_lo_se <- ldgr_se
brho_lapl_hi_lo_mean <- ldgr_mean

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
N_equil
N_brho_lo_hi <- 2227.849


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
N_pler_lo_hi <- 133.3819

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
N_lapl_lo_hi <- 116.9672


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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
pler_brho_lo_hi_se <- ldgr_se
pler_brho_lo_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
pler_lapl_lo_hi_se <- ldgr_se
pler_lapl_lo_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
lapl_brho_lo_hi_se <- ldgr_se
lapl_brho_lo_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
lapl_pler_lo_hi_se <- ldgr_se
lapl_pler_lo_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
brho_pler_lo_hi_se <- ldgr_se
brho_pler_lo_hi_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
brho_lapl_lo_hi_se <- ldgr_se
brho_lapl_lo_hi_mean <- ldgr_mean

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
N_equil
N_brho_lo_int <- 43.38967


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
N_pler_lo_int <- 164.6629

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
N_lapl_lo_int <- 15.91799


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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
pler_brho_lo_int_se <- ldgr_se
pler_brho_lo_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
pler_lapl_lo_int_se <- ldgr_se
pler_lapl_lo_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
lapl_brho_lo_int_se <- ldgr_se
lapl_brho_lo_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
lapl_pler_lo_int_se <- ldgr_se
lapl_pler_lo_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
brho_pler_lo_int_se <- ldgr_se
brho_pler_lo_int_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
brho_lapl_lo_int_se <- ldgr_se
brho_lapl_lo_int_mean <- ldgr_mean

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
N_pler_lo_lo <- 36.59878

##################
#######LAPL#######
##################
N_lapl_lo_lo <- 0

###########################
######CIs HI H2O INT N######
###########################
N_brho_lo_lo <- 0
N_lapl_lo_lo <- 0
N_pler_lo_lo <- 30.27704

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
pler_brho_lo_lo_se <- ldgr_se
pler_brho_lo_lo_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
pler_lapl_lo_lo_se <- ldgr_se
pler_lapl_lo_lo_mean <- ldgr_mean


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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
brho_pler_lo_lo_se <- ldgr_se
brho_pler_lo_lo_mean <- ldgr_mean

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
ldgr_se <- sd(ldgr_out)
ldgr_mean <- mean (ldgr_out)
brho_lapl_lo_lo_se <- ldgr_se
brho_lapl_lo_lo_mean <- ldgr_mean

grwr_mean_sd <- data.frame(brho_lapl_hi_hi_se,brho_lapl_hi_int_se,brho_lapl_hi_lo_se,
                      brho_lapl_lo_hi_se,brho_lapl_lo_int_se,brho_lapl_lo_lo_se,
                      pler_brho_hi_hi_se,pler_brho_hi_int_se,pler_brho_hi_lo_se,
                      pler_brho_lo_hi_se,pler_brho_lo_int_se,pler_brho_lo_lo_se,
                      lapl_brho_hi_hi_se,lapl_brho_hi_int_se,lapl_brho_lo_hi_se,
                      lapl_brho_lo_int_se,brho_pler_hi_hi_se,brho_pler_hi_int_se,
                      brho_pler_hi_lo_se,brho_pler_lo_hi_se,brho_pler_lo_int_se,
                      brho_pler_lo_lo_se,pler_lapl_hi_hi_se,pler_lapl_hi_int_se,
                      pler_lapl_hi_lo_se,pler_lapl_lo_hi_se,pler_lapl_lo_int_se,
                      pler_lapl_lo_lo_se,lapl_pler_hi_hi_se,lapl_pler_hi_int_se,
                      lapl_pler_lo_hi_se,lapl_pler_lo_int_se,brho_lapl_hi_hi_mean,brho_lapl_hi_int_mean,brho_lapl_hi_lo_mean,
                      brho_lapl_lo_hi_mean,brho_lapl_lo_int_mean,brho_lapl_lo_lo_mean,
                      pler_brho_hi_hi_mean,pler_brho_hi_int_mean,pler_brho_hi_lo_mean,
                      pler_brho_lo_hi_mean,pler_brho_lo_int_mean,pler_brho_lo_lo_mean,
                      lapl_brho_hi_hi_mean,lapl_brho_hi_int_mean,lapl_brho_lo_hi_mean,
                      lapl_brho_lo_int_mean,brho_pler_hi_hi_mean,brho_pler_hi_int_mean,
                      brho_pler_hi_lo_mean,brho_pler_lo_hi_mean,brho_pler_lo_int_mean,
                      brho_pler_lo_lo_mean,pler_lapl_hi_hi_mean,pler_lapl_hi_int_mean,
                      pler_lapl_hi_lo_mean,pler_lapl_lo_hi_mean,pler_lapl_lo_int_mean,
                      pler_lapl_lo_lo_mean,lapl_pler_hi_hi_mean,lapl_pler_hi_int_mean,
                      lapl_pler_lo_hi_mean,lapl_pler_lo_int_mean)

grwr_mean_sd2 <-pivot_longer(grwr_mean_sd,cols=1:64,names_to = "treatment", values_to = "sd") %>%
  mutate(ID = 1:64)

grwr_sd <- grwr_mean_sd2[grwr_mean_sd2$ID < 33, ] %>%
  select(1:2)

grwr_mean <- grwr_mean_sd2[grwr_mean_sd2$ID > 32, ] %>%
  rename("mean" = "sd") %>%
  select(1:2)

grwr_mean_sd_final <- cbind(grwr_sd,grwr_mean) %>%
  select(2:4)

write.csv(grwr_mean_sd_final,"grwr_sd_mean.csv")


