# example script for calculating LDGR with the posteriors

library(stan)
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

# extract posteriors
load(here("Stan models","no_dist_seeds_pler_hi_hi.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_hi.Rdata"))
pler_hi_hi <- rstan::extract(no_dist_seeds_pler_hi_hi)
brho_hi_hi <- rstan::extract(no_dist_seeds_brho_hi_hi)

##########################
###BRHO to equilibrium#####
##########################

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

##########################
###PLER invading BRHO#####
##########################
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
