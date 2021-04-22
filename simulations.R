library(tidyverse)

## Read in data
params2 <- read.csv(paste(datpath, "params2.csv", sep = ""))
trt <- read.csv(paste(datpath, "years_trt.csv", sep = "")) 
cover <-read.csv(paste(datpath,"JR_cover_1mplot.csv",sep=""))

## Determine equilibrium conditions for each species in isolation
params_dat <- params2

pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0)
  return(N)
}

## Bromus
bromus <-params_dat %>%
  filter(species == "Bromus")

N0 <- 27
time <- length(precipt$type_year)
N_bromus <- rep(NA, time)
N_bromus[1] <- N0
bs <- 0.013
bg <- 0.98

for (t in 1:time) {
  params <- subset(bromus,w_trt==precipt$type_year[t] & n_trt==nitrogen$n_trt[t])
  N_bromus[t+1] <- pop.equilibrium(N0=N_bromus[t], s=bs, g=bg, a_intra=bromus$alpha_Bromus, lambda=bromus$lambda)
}

# check output
plot(seq(1:(time+1)), N_bromus, type="l")
N_bromus

## Layia
layia <-params_dat %>%
  filter(species == "Layia")

N0 <- 89
time <- length(precipt$type_year)
N_layia <- rep(NA, time)
N_layia[1] <- N0
ls <- 0.15
lg <- 0.32

for (t in 1:time) {
  params <- subset(layia,w_trt==precipt$type_year[t] & n_trt==nitrogen$n_trt[t])
  N_layia[t+1] <- pop.equilibrium(N0=N_layia[t], s=ls, g=lg, a_intra=layia$alpha_Layia, lambda=layia$lambda)
}

# check output
plot(seq(1:(time+1)), N_layia, type="l")
N_layia

## Festuca
festuca <-params_dat %>%
  filter(species == "Festuca")

N0 <- 54
time <- length(precipt$type_year)
N_festuca <- rep(NA, time)
N_festuca[1] <- N0
fs <- 0.013
fg <- 0.83

for (t in 1:time) {
  params <- subset(festuca,w_trt==precipt$type_year[t] & n_trt==nitrogen$n_trt[t])
  N_festuca[t+1] <- pop.equilibrium(N0=N_festuca[t], s=fs, g=fg, a_intra=festuca$alpha_Festuca, lambda=festuca$lambda)
}

# check output
plot(seq(1:(time+1)), N_festuca, type="l")
N_festuca


## Data manipulation for starting conditions
# Equilibrium abundances
brho_eq <- 215
pler_eq <- 157
lapl_eq <- 152
femi_eq <- 485

params_dat <-params2 %>%
  separate(treatments,c("w_trt","n_trt"),sep="water.")

params_dat$w_trt[params_dat$w_trt == "hi."] <- "Wet"
params_dat$w_trt[params_dat$w_trt == "lo."] <- "Dry"

params_dat <- params_dat %>%
  unite(treatments,c(w_trt,n_trt),sep = ".")

trt <- trt %>%
  unite(treatments,c(type_year,n_trt),sep=".",remove=FALSE)
trt$n_trt[trt$n_trt == "lo.N"] <- "int.N"

cover_dat <- cover %>%
  filter(species == "BRMO" | species == "PLER" | species == "LAPL" | species == "VUMI") %>%
  group_by(year, species) %>%
  summarize(mean_cover = mean(cover)) %>%
  filter(year == 1983) %>%
  mutate(species = c("Bromus","Layia","Plantago","Festuca")) %>%
  mutate(abundance = mean_cover*(c(brho_eq,lapl_eq,pler_eq,femi_eq)/100))%>%
  select(-mean_cover) %>%
  pivot_wider(names_from = species, values_from = abundance)

cov_trt <- left_join(cover_dat,trt)
cov_trt_params <- left_join(cov_trt,params_dat)
trt_params <- left_join(trt,params_dat)

start_dat <- full_join(cov_trt_params,trt_params) 

#####################################
########## Simulation ###############
#####################################

start_dat <- start_dat %>%
  filter(year !=1983) %>%
  select(-Bromus,-Layia,-Plantago,-Festuca)

bg <- 0.98
bs <- 0.013
pg <- 0.92
ps <- 0.75
lg <- 0.32
ls <- 0.15
fg <- 0.83
fs <- 0.013

year <- length(1983:2019)

N = as.data.frame(matrix(NA, nrow=year, ncol=4))
colnames(N) = c("Nb", "Np", "Nf","Nl")
N[1,] =c(7.838542, 38.18153,27.31493,1.414444)


growth = function(N, start_dat,year){
  for (i in (year-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(start_dat$blambda[i]/(1 + start_dat$bap[i]*pg*N$Np[i] + start_dat$bab*bg*N$Nb[i] + start_dat$baf[i]*fg*N$Nf[i] + start_dat$bal[i]*lg*N$Nl[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pab[i]*bg*N$Nb[i]+ start_dat$paf[i]*fg*N$Nf[i] + start_dat$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$laf[i]*fg*N$Nf[i] + start_dat$lab[i]*bg*N$Nb[i] + start_dat$lap[i]*pg*N$Np[i]))
    
    N$Nf[i+1] = fs*(1-fg)*N$Nf[i] + fg*N$Nf[i]*(start_dat$flambda[i]/(1 + start_dat$faf[i]*fg*N$Nf[i] + start_dat$fal[i]*lg*N$Nl[i] + start_dat$fab[i]*bg*N$Nb[i] + start_dat$fap[i]*pg*N$Np[i]))
    
  }
  return(N)
}

growth(N,start_dat,year)
