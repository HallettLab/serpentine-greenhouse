library(tidyverse)
library(ggtext)

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
#femi_eq <- 485


trt <- read.csv(paste(datpath, "years_trt.csv", sep = "")) 

params_dat <-params2 %>%
  separate(treatments,c("w_trt","n_trt"),sep="water.")

params_dat$w_trt[params_dat$w_trt == "hi."] <- "Wet"
params_dat$w_trt[params_dat$w_trt == "lo."] <- "Dry"

params_dat <- params_dat %>%
  unite(treatments,c(w_trt,n_trt),sep = ".")

trt <- trt %>%
  unite(treatments,c(type_year,n_trt),sep=".",remove=FALSE)

cover_dat <- cover %>%
  filter(species == "BRMO" | species == "PLER" | species == "LAPL") %>%
  group_by(year, species) %>%
  summarize(mean_cover = mean(cover)) %>%
  filter(year == 1983) %>%
  mutate(species = c("Bromus","Layia","Plantago")) %>%
  mutate(abundance = mean_cover*c(brho_eq,lapl_eq,pler_eq)/100)%>%
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
  select(-Bromus,-Layia,-Plantago)

#write.csv(start_dat,"start_dat.csv")

bg <- 0.98
bs <- 0.013
pg <- 0.92
ps <- 0.75
lg <- 0.32
ls <- 0.15
#fg <- 0.83
#fs <- 0.013

year <- length(1983:2019)

N = as.data.frame(matrix(NA, nrow=37, ncol=4))
colnames(N) = c("Nb", "Np", "Nf","Nl")
N[1,] =c(7.838542, 38.18153,27.31493,1.414444)

# Nf: 27.31493, half: 13.65747, quarter: 6.828733

growth = function(N, start_dat,year){
  for (i in 1:(year-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(start_dat$blambda[i]/(1 + start_dat$bap[i]*pg*N$Np[i] + start_dat$bab[i]*bg*N$Nb[i] + start_dat$baf[i]*fg*N$Nf[i] + start_dat$bal[i]*lg*N$Nl[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pab[i]*bg*N$Nb[i]+ start_dat$paf[i]*fg*N$Nf[i] + start_dat$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$laf[i]*fg*N$Nf[i] + start_dat$lab[i]*bg*N$Nb[i] + start_dat$lap[i]*pg*N$Np[i]))
   
     N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
    N$Nf[i+1] = fs*(1-fg)*N$Nf[i] + fg*N$Nf[i]*(start_dat$flambda[i]/(1 + start_dat$faf[i]*fg*N$Nf[i] + start_dat$fal[i]*lg*N$Nl[i] + start_dat$fab[i]*bg*N$Nb[i] + start_dat$fap[i]*pg*N$Np[i]))
    
  }

  return(N)
}

N <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")
N$species[N$species == "Nb"] <- "Bromus"
N$species[N$species == "Np"] <- "Plantago"
N$species[N$species == "Nf"] <- "Festuca"
N$species[N$species == "Nl"] <- "Layia"

ggplot(N,aes(year,abundance,color=species)) + geom_line()

## Pairwise simulations
#BRHO and PLER
N = as.data.frame(matrix(NA, nrow=37, ncol=2))
colnames(N) = c("Nb", "Np")
N[1,] =c(7.838542, 38.18153)


growth = function(N, start_dat,year){
  for (i in 1:(year-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(start_dat$blambda[i]/(1 + start_dat$bap[i]*pg*N$Np[i] + start_dat$bab[i]*bg*N$Nb[i]))
    N$Nb[i+1] = ifelse(N$Nb[i+1]<0.1,0.1,N$Nb[i+1])
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pab[i]*bg*N$Nb[i]))
    
  }
  
  return(N)
}

N <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance") %>%
  mutate(cover = abundance/c(brho_eq/100,pler_eq/100))
N$species[N$species == "Nb"] <- "Bromus"
N$species[N$species == "Np"] <- "Plantago"

bp <- ggplot(N,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3) + xlab("Year") +
  theme(legend.text = element_markdown()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Plantago*"),
                     values=c("grey80","grey30"))

#PLER and LAPL
N = as.data.frame(matrix(NA, nrow=37, ncol=2))
colnames(N) = c("Np", "Nl")
N[1,] =c(38.18153,1.414444)

#Nl: 1.414444

growth = function(N, start_dat,year){
  for (i in 1:(year-1)){
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$lap[i]*pg*N$Np[i]))
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
  }
  
  return(N)
}

N <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")%>%
  mutate(cover = abundance/c(lapl_eq/100,pler_eq/100))
N$species[N$species == "Np"] <- "Plantago"
N$species[N$species == "Nl"] <- "Layia"

lp <- ggplot(N,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3) + xlab("Year") +
  theme(legend.text = element_markdown()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Layia*","*Plantago*"),
                     values=c("grey50","grey30"))

#BRHO,LAPL,PLER
N = as.data.frame(matrix(NA, nrow=37, ncol=3))
colnames(N) = c("Nb", "Np","Nl")
N[1,] =c(7.838542, 38.18153,1.414444)

#Nl:1.414444

growth = function(N, start_dat,year){
  for (i in 1:(year-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(start_dat$blambda[i]/(1 + start_dat$bap[i]*pg*N$Np[i] + start_dat$bab[i]*bg*N$Nb[i] + start_dat$bal[i]*lg*N$Nl[i]))
    N$Nb[i+1] = ifelse(N$Nb[i+1]<0.1,0.1,N$Nb[i+1])
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pab[i]*bg*N$Nb[i]+ start_dat$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$lab[i]*bg*N$Nb[i] + start_dat$lap[i]*pg*N$Np[i]))
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
    
  }
  
  return(N)
}

N <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")%>%
  mutate(cover = abundance/c(brho_eq/100,lapl_eq/100,pler_eq/100))
N$species[N$species == "Nb"] <- "Bromus"
N$species[N$species == "Np"] <- "Plantago"
N$species[N$species == "Nl"] <- "Layia"

blp <- ggplot(N,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3) + xlab("Year") +
  theme(legend.text = element_markdown()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("grey80","grey50","grey30"))
#BRHO, PLER, and LAPL
brho_bg <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98)
start_dat_brho <- start_dat %>%
  cbind(bg = brho_bg)

N = as.data.frame(matrix(NA, nrow=37, ncol=3))
colnames(N) = c("Nb", "Np","Nl")
N[1,] =c(7.838542, 38.18153,1.414444)

#Nl:1.414444

growth = function(N, start_dat,year){
  for (i in 1:(year-1)){
    N$Nb[i+1] = bs*(1-start_dat_brho$bg[i])*N$Nb[i]  + start_dat_brho$bg[i]*N$Nb[i]*(start_dat$blambda[i]/(1 + start_dat$bap[i]*pg*N$Np[i] + start_dat$bab[i]*start_dat_brho$bg[i]*N$Nb[i] + start_dat$bal[i]*lg*N$Nl[i]))
    N$Nb[i+1] = ifelse(N$Nb[i+1]<0.1,0.1,N$Nb[i+1])
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pab[i]*start_dat_brho$bg[i]*N$Nb[i]+ start_dat$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$lab[i]*start_dat_brho$bg[i]*N$Nb[i] + start_dat$lap[i]*pg*N$Np[i]))
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    

  }
  
  return(N)
}

N <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")%>%
  mutate(cover = abundance/c(brho_eq/100,lapl_eq/100,pler_eq/100))
N$species[N$species == "Nb"] <- "Bromus"
N$species[N$species == "Np"] <- "Plantago"
N$species[N$species == "Nl"] <- "Layia"

blp2 <- ggplot(N,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3) + xlab("Year") +
  theme(legend.text = element_markdown()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("grey80","grey50","grey30"))

# Bromus and Layia
N = as.data.frame(matrix(NA, nrow=37, ncol=2))
colnames(N) = c("Nb","Nl")
N[1,] =c(7.838542,1.414444)

#Nl:1.414444

growth = function(N, start_dat,year){
  for (i in 1:(year-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(start_dat$blambda[i]/(1 + start_dat$bab[i]*bg*N$Nb[i] + start_dat$bal[i]*lg*N$Nl[i]))
    N$Nb[i+1] = ifelse(N$Nb[i+1]<0.1,0.1,N$Nb[i+1])
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$lab[i]*bg*N$Nb[i]))
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
    
  }
  
  return(N)
}

N <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")
N$species[N$species == "Nb"] <- "Bromus"
N$species[N$species == "Nl"] <- "Layia"

bl <- ggplot(N,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3) + xlab("Year") +
  theme(legend.text = element_markdown()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*"),
                     values=c("grey80","grey50"))

#arrange plots
ggarrange(lp,jr,bp,blp,bl,blp2,
          labels = c("a)", "d)", "b)","e)","c)","f)"),
          ncol = 2, nrow = 3)

