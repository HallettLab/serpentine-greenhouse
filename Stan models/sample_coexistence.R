
## Full annual plant population models
### Run single species to equilibrium conexistence
run.to.equilibrium <- function(surv, germ, lambda, alpha_intra, Nt) {
  Ntp1 <- (1-germ)*surv*Nt + germ*lambda*Nt/(1+ alpha_intra * Nt)
  return(Ntp1)
}

### Run invader into resident at equilibrium
run.invader <- function(surv, germ, lambda, alpha_inter, resid_abund, invader_abund) {
  Ntp1 <- (1-germ)*surv*invader_abund + germ*lambda*invader_abund/(1+ alpha_inter * resid_abund)
  LDGR <- log(Ntp1/invader_abund)
  return(LDGR)
}

## Set germination and survival rates estimated from lit
# Gulmon (germ) and LaForgia (surv)
pler_hi_hi$germ <- .92
pler_hi_int$germ <- .92
pler_hi_lo$germ <- .92
pler_lo_hi$germ <- .92
pler_lo_int$germ <- .92
pler_lo_lo$germ <- .92
pler_hi_hi$surv <- .75
pler_hi_int$surv <- .75
pler_hi_lo$surv <- .75
pler_lo_hi$surv <- .75
pler_lo_int$surv <- .75
pler_lo_lo$surv <- .75

# Gulmon (germ) and Andrew's code (surv)
brho_hi_hi$germ <- .98
brho_hi_int$germ <- .98
brho_hi_lo$germ <- .98
brho_lo_hi$germ <- .98
brho_lo_int$germ <- .98
brho_lo_lo$germ <- .98
brho_hi_hi$surv <- .013
brho_hi_int$surv <- .013
brho_hi_lo$surv <- .013
brho_lo_hi$surv <- .013
brho_lo_int$surv <- .013
brho_lo_lo$surv <- .013

# Rossington M.S. thesis (germ and surv)
lapl_hi_hi$germ <- .32
lapl_hi_int$germ <- .32
lapl_hi_lo$germ <- .32
lapl_lo_hi$germ <- .32
lapl_lo_int$germ <- .32
lapl_lo_lo$germ <- .32
lapl_hi_hi$surv <- .32
lapl_hi_int$surv <- .15
lapl_hi_lo$surv <- .15
lapl_lo_hi$surv <- .15
lapl_lo_int$surv <- .15
#lapl_lo_lo$surv <- .15

# Gulmon (germ) and X (surv)
#FEMI surv based on BRHO surv
femi_hi_hi$germ <- .83
femi_hi_int$germ <- .83
femi_hi_lo$germ <- .83
femi_lo_hi$germ <- .83
femi_lo_int$germ <- .83
femi_lo_lo$germ <- .83
femi_hi_hi$surv <- .013
femi_hi_int$surv <- .013
femi_hi_lo$surv <- .013
femi_lo_hi$surv <- .013
femi_lo_int$surv <- .013
femi_lo_lo$surv <- .013

## Run to equilibrium each species alone in each treatment
all_datset <- list(pler_hi_hi,pler_hi_int,pler_hi_lo,pler_lo_hi,pler_lo_int,pler_lo_lo,
                   brho_hi_hi, brho_hi_int, brho_hi_lo, brho_lo_hi,brho_lo_int,brho_lo_lo) # set which dataset we want to do
all_intra <- c("alpha_pler", "alpha_pler", "alpha_brho", "alpha_brho")
options <- length(all_intra)

time <- 200
runs <- 200

N <- array(NA, c(options, runs, time))
N[,,1] <- 100

### Loop once for each species x treatment 
for (x in 1:options) {
  datset <- all_datset[[x]]
  intra <- all_intra[[x]]
  
  post_length <- length(datset$lambda)
  
  all_intras <- datset[[intra]]
  
  for (t in 1:(time-1)) {
    posts <- sample(post_length, runs, replace=TRUE)
    lambda <- datset$lambda[posts]
    alpha_intra <- all_intras[posts]
    N[x, ,t+1] <- run.to.equilibrium(surv=datset$surv, germ=datset$germ, 
                                     lambda=lambda, alpha_intra=alpha_intra, Nt=N[x, ,t])
  }
}

## Split into two objects, one for wet and one for dry treatments, and add labels
residents_hi_hi <-  data.frame(N[1,,200], N[3,,200], N[5,,200], N[7,,200], N[9,,200])
names(residents_dry) <- c("PLER", "BRHO")

residents_hi_int <-  data.frame(N[1,,200], N[3,,200], N[5,,200], N[7,,200], N[9,,200])
names(residents_wet) <- c("PLER", "BRHO")

residents_hi_lo <-  data.frame(N[1,,200], N[3,,200], N[5,,200], N[7,,200], N[9,,200])
names(residents_dry) <- c("PLER", "BRHO")

residents_lo_hi <-  data.frame(N[1,,200], N[3,,200], N[5,,200], N[7,,200], N[9,,200])
names(residents_wet) <- c("PLER", "BRHO")

residents_lo_int <-  data.frame(N[1,,200], N[3,,200], N[5,,200], N[7,,200], N[9,,200])
names(residents_dry) <- c("PLER", "BRHO")

residents_lo_lo <-  data.frame(N[2,,200], N[4,,200], N[6,,200], N[8,,200], N[10,,200])
names(residents_wet) <- c("PLER", "BRHO")



## Set up objects to store results of Avena invasions into each other species in the dry treatment 
reps <- 200
avfa_into_brho_dry <- matrix(NA, reps, runs)
avfa_into_vumy_dry <- matrix(NA, reps, runs)
avfa_into_laca_dry <- matrix(NA, reps, runs)
avfa_into_esca_dry <- matrix(NA, reps, runs)

## Invade Avena into each species in the dry treatment

invader_abund <- 1

post_length <- length(avfa_dry$lambda)
for (r in 1:reps) {
  
  posts <- sample(post_length, runs, replace=TRUE)
  
  avfa_into_brho_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_brho[posts],
                                        resid_abund=residents_dry$BRHO, invader_abund=invader_abund)
  
  avfa_into_vumy_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_vumy[posts],
                                        resid_abund=residents_dry$VUMY, invader_abund=invader_abund)
  
  avfa_into_laca_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_laca[posts],
                                        resid_abund=residents_dry$LACA, invader_abund=invader_abund)
  
  avfa_into_esca_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_esca[posts],
                                        resid_abund=residents_dry$ESCA, invader_abund=invader_abund)
}
