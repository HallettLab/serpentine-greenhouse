library(tidyverse)

## Read in data
datpath <- "C:/Users/ehernan2/Dropbox (University Of Oregon)/Thesis/"
params <- read.csv(paste(datpath, "params.csv", sep = ""))
trt <- read.csv(paste(datpath, "years_trt.csv", sep = "")) 
cover <-read.csv(paste(datpath,"JR_cover_1mplot.csv",sep=""))

## Determine equilibrium conditions for each species in isolation
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

## Equilibrium abundances
brho_eq <- 215
pler_eq <- 157
lapl_eq <- 152
femi_eq <- 485


## Data manipulation for starting conditions
params_dat <-params

params_dat$germ[params_dat$species == "Bromus"] <- 0.98
params_dat$germ[params_dat$species == "Festuca"] <- 0.83
params_dat$germ[params_dat$species == "Plantago"] <- 0.92
params_dat$germ[params_dat$species == "Layia"] <- 0.32
params_dat$surv[params_dat$species == "Bromus"] <- 0.013
params_dat$surv[params_dat$species == "Festuca"] <- 0.013
params_dat$surv[params_dat$species == "Plantago"] <- 0.75
params_dat$surv[params_dat$species == "Layia"] <- 0.15
names(params_dat)[3] <- "alpha_Bromus"
names(params_dat)[4] <- "alpha_Festuca"
names(params_dat)[5] <- "alpha_Layia"
names(params_dat)[6] <- "alpha_Plantago"

params_dat <- params_dat %>%
  separate(treatments,c("w_trt","n_trt"),sep="water.") %>%
  filter(species != "Bromus2")

params_dat$w_trt[params_dat$w_trt == "hi."] <- "Wet"
params_dat$w_trt[params_dat$w_trt == "lo."] <- "Dry"

params_dat <- params_dat %>%
  unite(treatments,c(w_trt,n_trt),sep = ".")

trt2 <- trt %>%
  unite(treatments,c("type_year","n_trt"),sep=".",remove=FALSE)%>%
  filter(year!=1983) %>%
  slice(rep(1:n(), each = 4))

species<-as.data.frame(rep(c("Bromus","Layia","Plantago","Festuca"),times=36))
colnames(species) = "species"

trt3 <-trt2 %>%
cbind(species)

cover_dat <- cover %>%
  filter(species == "BRMO" | species == "PLER" | species == "LAPL" | species == "VUMI") %>%
  group_by(year, species) %>%
  summarize(mean_cover = mean(cover)) %>%
  filter(year == 1983) %>%
  mutate(species = c("Bromus","Layia","Plantago","Festuca"))

cov_trt <- left_join(cover_dat,trt)
cov_trt_params <- left_join(cov_trt,params_dat)
trt_params <- left_join(trt3,params_dat)

start_dat <- full_join(cov_trt_params,trt_params) %>%
  mutate(abundance = mean_cover*(c(brho_eq,lapl_eq,pler_eq,femi_eq)/100))
start_dat <- start_dat[,c(1,2,3,15,4,5,6,7,8,9,10,11,12,13,14)]

write.csv(start_dat,"start_sim.csv")

