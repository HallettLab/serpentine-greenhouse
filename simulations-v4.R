library(tidyverse)
library(ggtext)
library(cowplot)
library(grid)
library(ggplotify)
library(rstan)
library(here)

ps <- .75 # gulmon
pg <- .92 # gulmon
bs <- .01 # andrew
bg <- .98 # gulmon
ls <- .15 # rossington
lg <- .32 # rossington

## Read in data
#params2 <- read.csv(paste(datpath, "params2_LGS.csv", sep = "")) # parameters from first stan model fits
trt <- read.csv(paste(datpath, "years_trt.csv", sep = "")) %>%
  rename(w_trt=type_year) %>%
  filter(year %in% 1983:2019)
cover <-read.csv(paste(datpath,"JR_cover_1mplot1983-2019.csv",sep=""))
#stems_dat <-  read.csv(paste(datpath, "/stems_background.csv", sep = "")) %>%
  #filter(seed_sp != "FEMI") %>%
  #mutate(recruit=stem_density/seed_added) 

##survival and germination fractions
ps <- .75 # gulmon
pg <- .92 # gulmon
bs <- .01 # andrew
bg <- .98 # gulmon
ls <- .15 # rossington
lg <- .32 # rossington

##Load stan models from Stan models folder
load(here("Stan models_3spp","pler_hi_hi.Rdata"))
load(here("Stan models_3spp","brho_hi_hi.Rdata"))
load(here("Stan models_3spp","lapl_hi_hi.Rdata"))
load(here("Stan models_3spp","pler_hi_int.Rdata"))
load(here("Stan models_3spp","brho_hi_int.Rdata"))
load(here("Stan models_3spp","lapl_hi_int.Rdata"))
load(here("Stan models_3spp","pler_hi_lo.Rdata"))
load(here("Stan models_3spp","brho_hi_lo.Rdata"))
load(here("Stan models_3spp","pler_lo_hi.Rdata"))
load(here("Stan models_3spp","brho_lo_hi.Rdata"))
load(here("Stan models_3spp","lapl_lo_hi.Rdata"))
load(here("Stan models_3spp","pler_lo_int.Rdata"))
load(here("Stan models_3spp","brho_lo_int.Rdata"))
load(here("Stan models_3spp","lapl_lo_int.Rdata"))
load(here("Stan models_3spp","pler_lo_lo.Rdata"))
load(here("Stan models_3spp","brho_lo_lo.Rdata"))

#Extract parameters from models and rename 
brho_hi_hi <- rstan::extract(brho_hi_hi)
brho_hi_int <- rstan::extract(brho_hi_int)
brho_hi_lo <- rstan::extract(brho_hi_lo)
brho_lo_hi <- rstan::extract(brho_lo_hi)
brho_lo_int <- rstan::extract(brho_lo_int)
brho_lo_lo <- rstan::extract(brho_lo_lo)
lapl_hi_hi <- rstan::extract(lapl_hi_hi)
lapl_hi_int <- rstan::extract(lapl_hi_int)
lapl_lo_hi <- rstan::extract(lapl_lo_hi)
lapl_lo_int <- rstan::extract(lapl_lo_int)
pler_hi_hi <- rstan::extract(pler_hi_hi)
pler_hi_int <- rstan::extract(pler_hi_int)
pler_hi_lo <- rstan::extract(pler_hi_lo)
pler_lo_hi <- rstan::extract(pler_lo_hi)
pler_lo_int <- rstan::extract(pler_lo_int)
pler_lo_lo <- rstan::extract(pler_lo_lo)

##hi water hi N posteriors
# brho
posterior_length <- length(brho_hi_hi$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
blambda <- brho_hi_hi$lambda[posts]
bab <- brho_hi_hi$alpha_brho[posts]
bap <- brho_hi_hi$alpha_pler[posts]
bal <- brho_hi_hi$alpha_lapl[posts]
#lapl
posterior_length <- length(lapl_hi_hi$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
llambda <- lapl_hi_hi$lambda[posts]
lal <- lapl_hi_hi$alpha_lapl[posts]
lap <- lapl_hi_hi$alpha_pler[posts]
lab <- lapl_hi_hi$alpha_brho[posts]
#pler
posterior_length <- length(pler_hi_hi$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
plambda <- pler_hi_hi$lambda[posts]
pap <- pler_hi_hi$alpha_pler[posts]
pab <- pler_hi_hi$alpha_brho[posts]
pal <- pler_hi_hi$alpha_lapl[posts]

hi_hi_posts <-as.data.frame(cbind(blambda,bab,bap,bal,llambda,lal,lap,lab,plambda,pap,pab,pal))%>%
  mutate(w_trt="Wet",n_trt="hi.N")

##hi water int N posteriors
# brho
posterior_length <- length(brho_hi_int$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
blambda <- brho_hi_int$lambda[posts]
bab <- brho_hi_int$alpha_brho[posts]
bap <- brho_hi_int$alpha_pler[posts]
bal <- brho_hi_int$alpha_lapl[posts]
#lapl
posterior_length <- length(lapl_hi_int$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
llambda <- lapl_hi_int$lambda[posts]
lal <- lapl_hi_int$alpha_lapl[posts]
lap <- lapl_hi_int$alpha_pler[posts]
lab <- lapl_hi_int$alpha_brho[posts]
#pler
posterior_length <- length(pler_hi_int$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
plambda <- pler_hi_int$lambda[posts]
pap <- pler_hi_int$alpha_pler[posts]
pab <- pler_hi_int$alpha_brho[posts]
pal <- pler_hi_int$alpha_lapl[posts]

hi_int_posts <-as.data.frame(cbind(blambda,bab,bap,bal,llambda,lal,lap,lab,plambda,pap,pab,pal))%>%
  mutate(w_trt="Wet",n_trt="int.N")

##hi water lo N posteriors
# brho
posterior_length <- length(brho_hi_lo$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
blambda <- brho_hi_lo$lambda[posts]
bab <- brho_hi_lo$alpha_brho[posts]
bap <- brho_hi_lo$alpha_pler[posts]
bal <- brho_hi_lo$alpha_lapl[posts]
#lapl
llambda <- 0
lal <- 0
lap <- 0
lab <- 0
#pler
posterior_length <- length(pler_hi_lo$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
plambda <- pler_hi_lo$lambda[posts]
pap <- pler_hi_lo$alpha_pler[posts]
pab <- pler_hi_lo$alpha_brho[posts]
pal <- pler_hi_lo$alpha_lapl[posts]

hi_lo_posts <-as.data.frame(cbind(blambda,bab,bap,bal,llambda,lal,lap,lab,plambda,pap,pab,pal))%>%
  mutate(w_trt="Wet",n_trt="lo.N")

##lo water hi N posteriors
# brho
posterior_length <- length(brho_lo_hi$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
blambda <- brho_lo_hi$lambda[posts]
bab <- brho_lo_hi$alpha_brho[posts]
bap <- brho_lo_hi$alpha_pler[posts]
bal <- brho_lo_hi$alpha_lapl[posts]
#lapl
posterior_length <- length(lapl_lo_hi$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
llambda <- lapl_lo_hi$lambda[posts]
lal <- lapl_lo_hi$alpha_lapl[posts]
lap <- lapl_lo_hi$alpha_pler[posts]
lab <- lapl_lo_hi$alpha_brho[posts]
#pler
posterior_length <- length(pler_lo_hi$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
plambda <- pler_lo_hi$lambda[posts]
pap <- pler_lo_hi$alpha_pler[posts]
pab <- pler_lo_hi$alpha_brho[posts]
pal <- pler_lo_hi$alpha_lapl[posts]

lo_hi_posts <-as.data.frame(cbind(blambda,bab,bap,bal,llambda,lal,lap,lab,plambda,pap,pab,pal))%>%
  mutate(w_trt="Dry",n_trt="hi.N")

##lo water int N posteriors
# brho
posterior_length <- length(brho_lo_int$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
blambda <- brho_lo_int$lambda[posts]
bab <- brho_lo_int$alpha_brho[posts]
bap <- brho_lo_int$alpha_pler[posts]
bal <- brho_lo_int$alpha_lapl[posts]
#lapl
posterior_length <- length(lapl_lo_int$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
llambda <- lapl_lo_int$lambda[posts]
lal <- lapl_lo_int$alpha_lapl[posts]
lap <- lapl_lo_int$alpha_pler[posts]
lab <- lapl_lo_int$alpha_brho[posts]
#pler
posterior_length <- length(pler_lo_int$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
plambda <- pler_lo_int$lambda[posts]
pap <- pler_lo_int$alpha_pler[posts]
pab <- pler_lo_int$alpha_brho[posts]
pal <- pler_lo_int$alpha_lapl[posts]

lo_int_posts <-as.data.frame(cbind(blambda,bab,bap,bal,llambda,lal,lap,lab,plambda,pap,pab,pal))%>%
  mutate(w_trt="Dry",n_trt="int.N")

##lo water lo N posteriors
# brho
posterior_length <- length(brho_lo_lo$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
blambda <- brho_lo_lo$lambda[posts]
bab <- brho_lo_lo$alpha_brho[posts]
bap <- brho_lo_lo$alpha_pler[posts]
bal <- brho_lo_lo$alpha_lapl[posts]
#lapl
llambda <- 0
lal <- 0
lap <- 0
lab <- 0
#pler
posterior_length <- length(pler_lo_lo$lambda)
posts <- sample(posterior_length, 2000, replace=FALSE)
plambda <- pler_lo_lo$lambda[posts]
pap <- pler_lo_lo$alpha_pler[posts]
pab <- pler_lo_lo$alpha_brho[posts]
pal <- pler_lo_lo$alpha_lapl[posts]

lo_lo_posts <-as.data.frame(cbind(blambda,bab,bap,bal,llambda,lal,lap,lab,plambda,pap,pab,pal))%>%
  mutate(w_trt="Dry",n_trt="lo.N")

posts <- as.data.frame(rbind(hi_hi_posts,hi_int_posts,hi_lo_posts,lo_hi_posts,lo_int_posts,lo_lo_posts))

trt_posts <- full_join(trt,posts)

##join posts with trt dataframe
cover_dat <- cover %>%
  filter(treatment == "c") %>%
  filter(species == "PLER" | species == "BRMO" | species == "LAPL") %>%
  filter(year==1983) %>%
  group_by(species) %>%
  summarize(mean_cover=mean(cover)/100) %>%
  pivot_wider(names_from=species,values_from=mean_cover) %>%
  mutate(year=1983)

source("equil-abund-by-trt.R")

#1983 is wet and low N
species_eq_wet_lo <- species_eq %>%
  filter(w_trt=="Wet") %>%
  filter(n_trt=="lo.N") 

BRHO_equil <- species_eq_wet_lo %>%
  filter(species=="Bromus")

LAPL_equil <- species_eq_wet_lo %>%
  filter(species=="Layia")

PLER_equil <- species_eq_wet_lo %>%
  filter(species=="Plantago")

cover_trt <- left_join(cover_dat,trt) %>%
  rename(BRHO = BRMO) %>%
  mutate(BRHO = BRHO*BRHO_equil$equil_abundance,PLER=PLER*PLER_equil$equil_abundance,LAPL=LAPL*LAPL_equil$equil_abundance)

alldat0 <- left_join(trt_posts,cover_trt) %>%
  filter(year !=1983) %>%
  select(-BRHO, -LAPL, -PLER)

key <- as.data.frame(cbind(replicate.var=rep(seq(1:2000),36)))

alldat <- cbind(key, alldat0) %>%
  mutate(Nb = cover_trt$BRHO, Nl=cover_trt$LAPL, Np=cover_trt$PLER)

#dummy <- alldat %>%
#filter(replicate.var == 1)

######################################
########## Simulations ###############
######################################

###all species sim 1 without altering bromus germination###

growth = function(N){

  for (i in 1:(nrow(N)-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(N$blambda[i]/(1 + N$bap[i]*pg*N$Np[i] + N$bab[i]*bg*N$Nb[i] + N$bal[i]*lg*N$Nl[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(N$plambda[i]/(1 + N$pap[i]*pg*N$Np[i] + N$pab[i]*bg*N$Nb[i] + N$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(N$llambda[i]/(1 + N$lal[i]*lg*N$Nl[i] + N$lab[i]*bg*N$Nb[i] + N$lap[i]*pg*N$Np[i]))
   
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
    
  }

  return(N)
}

#growth(dummy)

X <- split(alldat, alldat["replicate.var"])
out <- lapply(X, FUN=growth)
output_sim1 <- do.call("rbind",out)

## data manipulation
sim1_mean <- output_sim1 %>%
  group_by(year) %>%
  summarize(Bromus=mean(Nb),Layia=mean(Nl),Plantago=mean(Np)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "mean_abundance")
sim1_CI25 <- output_sim1 %>%
  group_by(year) %>%
  summarize(Bromus=quantile(Nb,probs = 0.25),Layia=quantile(Nl,probs = 0.25),Plantago=quantile(Np,probs = 0.25)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "CIlower") 
sim1_CI75 <- output_sim1 %>%
  group_by(year) %>%
  summarize(Bromus=quantile(Nb,probs = 0.75),Layia=quantile(Nl,probs = 0.75),Plantago=quantile(Np,probs = 0.75)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "CIupper") 
sim1_CI <- left_join(sim1_CI25,sim1_CI75)
sim1_med <- output_sim1 %>%
  group_by(year) %>%
  summarize(Bromus=median(Nb),Layia=median(Nl),Plantago=median(Np)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "median") 
sim1 <- left_join(sim1_mean,sim1_CI) 
sim1 <- left_join(sim1,sim1_med)
cols <- c("mean_abundance","CIlower","CIupper","median")
sim1[cols] <- log(sim1[cols])
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
sim1[is.nan(sim1)] <- 0
sim1["mean_abundance"][sim1["mean_abundance"] == -Inf] <- 0
sim1["median"][sim1["median"] == -Inf] <- 0
sim1["CIlower"][sim1["CIlower"] == -Inf] <- 0
sim1["CIupper"][sim1["CIupper"] == -Inf] <- 0
sim1["median"][sim1["median"] == -Inf] <- 0
sim1$mean_abundance <- ifelse(sim1$mean_abundance < 0, 0, sim1$mean_abundance)
sim1$CIlower <- ifelse(sim1$CIlower < 0, 0, sim1$CIlower)
sim1$CIupper <- ifelse(sim1$CIupper < 0, 0, sim1$CIupper)
sim1$median <- ifelse(sim1$median < 0, 0, sim1$median)

########################################################
## all species sim 2 with changes in Bromus germination ##
########################################################
bg2 <- rep(0.98,72000)
alldatbg <- cbind(bg2,alldat)
colnames(alldatbg)[1] ="bg2"



growth = function(N){
  
  for (i in 1:(nrow(N)-1)){
    N$Nb[i+1] = bs*(1-N$bg2[i])*N$Nb[i]  + N$bg2[i]*N$Nb[i]*(N$blambda[i]/(1 + N$bap[i]*pg*N$Np[i] + N$bab[i]*N$bg2[i]*N$Nb[i] + N$bal[i]*lg*N$Nl[i]))
    
    N$Nb[1:26] = 0
    
    N$Nb[27] = 2
    
    N$Nb[28] = (bs*(1-N$bg2[27])*2  + N$bg2[27]*2*(N$blambda[27]/(1 + N$bap[27]*pg*N$Np[27] + N$bab[27]*N$bg2[27]*2 + N$bal[27]*lg*N$Nl[27]))) + 4
    
    N$Nb[29] = (bs*(1-N$bg2[28])*N$Nb[28]  + N$bg2[28]*N$Nb[28]*(N$blambda[28]/(1 + N$bap[28]*pg*N$Np[28] + N$bab[28]*N$bg2[28]*N$Nb[28] + N$bal[28]*lg*N$Nl[28]))) + 4
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(N$plambda[i]/(1 + N$pap[i]*pg*N$Np[i] + N$pab[i]*N$bg2[i]*N$Nb[i] + N$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(N$llambda[i]/(1 + N$lal[i]*lg*N$Nl[i] + N$lab[i]*N$bg2[i]*N$Nb[i] + N$lap[i]*pg*N$Np[i]))
    
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
    
  }
  
  return(N)
}

X <- split(alldatbg, alldat["replicate.var"])
out <- lapply(X, FUN=growth)
output_sim2 <- do.call("rbind",out)

sim2_mean <- output_sim2 %>%
  group_by(year) %>%
  summarize(Bromus=mean(Nb),Layia=mean(Nl),Plantago=mean(Np)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "mean_abundance")
sim2_CI25 <- output_sim2 %>%
  group_by(year) %>%
  summarize(Bromus=quantile(Nb,probs = 0.25),Layia=quantile(Nl,probs = 0.25),Plantago=quantile(Np,probs = 0.25)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "CIlower") 
sim2_CI75 <- output_sim2 %>%
  group_by(year) %>%
  summarize(Bromus=quantile(Nb,probs = 0.75),Layia=quantile(Nl,probs = 0.75),Plantago=quantile(Np,probs = 0.75)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "CIupper") 
sim2_CI <- left_join(sim2_CI25,sim2_CI75)
sim2_med <- output_sim2 %>%
  group_by(year) %>%
  summarize(Bromus=median(Nb),Layia=median(Nl),Plantago=median(Np)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "median") 
sim2 <- left_join(sim2_mean,sim2_CI) 
sim2 <- left_join(sim2,sim2_med)
cols <- c("mean_abundance","CIlower","CIupper","median")
sim2[cols] <- log(sim2[cols])
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
sim2[is.nan(sim2)] <- 0
sim2["mean_abundance"][sim2["mean_abundance"] == -Inf] <- 0
sim2["median"][sim2["median"] == -Inf] <- 0
sim2["CIlower"][sim2["CIlower"] == -Inf] <- 0
sim2["CIupper"][sim2["CIupper"] == -Inf] <- 0
sim2["median"][sim2["median"] == -Inf] <- 0
sim2$mean_abundance <- ifelse(sim2$mean_abundance < 0, 0, sim2$mean_abundance)
sim2$CIlower <- ifelse(sim2$CIlower < 0, 0, sim2$CIlower)
sim2$CIupper <- ifelse(sim2$CIupper < 0, 0, sim2$CIupper)
sim2$median <- ifelse(sim2$median < 0, 0, sim2$median)

ggplot(sim2,aes(year,median,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  geom_ribbon(aes(x= year, ymin=CIlower,ymax=CIupper, fill=species), color=NA, alpha = .2) +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"),legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_blank()) +
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#48B99B")) +
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                    values=c("#D55E00","#0072B2","#48B99B")) +
  scale_x_continuous(expand = c(0.04,0.04))


#############################
######PLER and LAPL sim######
#############################
growth = function(N){
  
  for (i in 1:(nrow(N)-1)){
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(N$plambda[i]/(1 + N$pap[i]*pg*N$Np[i] + N$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(N$llambda[i]/(1 + N$lal[i]*lg*N$Nl[i] + N$lap[i]*pg*N$Np[i]))
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
    
  }
  
  return(N)
}


X <- split(alldat, alldat["replicate.var"])
out <- lapply(X, FUN=growth)
output_simlp <- do.call("rbind",out)

simlp_mean <- output_simlp %>%
  group_by(year) %>%
  summarize(Bromus=mean(Nb),Layia=mean(Nl),Plantago=mean(Np)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "mean_abundance")
simlp_CI25 <- output_simlp %>%
  group_by(year) %>%
  summarize(Bromus=quantile(Nb,probs = 0.25),Layia=quantile(Nl,probs = 0.25),Plantago=quantile(Np,probs = 0.25)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "CIlower") 
simlp_CI75 <- output_simlp %>%
  group_by(year) %>%
  summarize(Bromus=quantile(Nb,probs = 0.75),Layia=quantile(Nl,probs = 0.75),Plantago=quantile(Np,probs = 0.75)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "CIupper") 
simlp_CI <- left_join(simlp_CI25,simlp_CI75)
simlp_med <- output_simlp %>%
  group_by(year) %>%
  summarize(Bromus=median(Nb),Layia=median(Nl),Plantago=median(Np)) %>%
  pivot_longer(2:4,names_to = "species",values_to = "median") 
simlp <- left_join(simlp_mean,simlp_CI) 
simlp <- left_join(simlp,simlp_med)
cols <- c("mean_abundance","CIlower","CIupper","median")
simlp[cols] <- log(simlp[cols])
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
simlp[is.nan(simlp)] <- 0
simlp["mean_abundance"][simlp["mean_abundance"] == -Inf] <- 0
simlp["median"][simlp["median"] == -Inf] <- 0
simlp["CIlower"][simlp["CIlower"] == -Inf] <- 0
simlp["CIupper"][simlp["CIupper"] == -Inf] <- 0
simlp["median"][simlp["median"] == -Inf] <- 0
simlp$mean_abundance <- ifelse(simlp$mean_abundance < 0, 0, simlp$mean_abundance)
simlp$CIlower <- ifelse(simlp$CIlower < 0, 0, simlp$CIlower)
simlp$CIupper <- ifelse(simlp$CIupper < 0, 0, simlp$CIupper)
simlp$median <- ifelse(simlp$median < 0, 0, simlp$median)

################################
#####Figure visualization#######
################################
theme_set(theme_bw())
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 16),
              axis.text = element_text(size = 12))

trts <- ggplot(trt,aes(year,growing_season_ppt)) +
  annotate("rect", xmin = -Inf, xmax = 1995, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow2") +
  annotate("rect", xmin = 1995, xmax = 2007, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow3")+
  annotate("rect", xmin = 2007, xmax = Inf, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow4")+
  geom_line(size=0.8)+
  geom_hline(yintercept=565,lty=2)+
  geom_point(size=3, aes(shape=w_trt))+
  scale_shape_manual(name="Growing season type",values=c(1,16)) +
  xlab("Year") +
  ylab("Precipitation (mm)")+
  theme(legend.position = "top",axis.text.x = element_blank(),axis.title.x = element_blank(),
        plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"))+
  scale_x_continuous(expand = c(0.04,0.04)) +
  annotate("text", x=1988,y=1200, label="Low N")+
  annotate("text", x=2001,y=1200, label="Intermediate N") +
  annotate("text", x=2014,y=1200, label="High N") 

trt2 <- trt %>%
  mutate(gst = 1,n=NA) 

trt2$n[trt2$n_trt == "lo.N"] <- "Low"
trt2$n[trt2$n_trt == "int.N"] <- "Intermediate"
trt2$n[trt2$n_trt == "hi.N"] <- "High"

trt2$n <- factor(trt2$n,levels = c("Low","Intermediate","High"))

bar <- ggplot(trt2,aes(year,gst)) +
  geom_point(aes(year,gst,color=n),size=0,shape =15)+
  annotate("rect", xmin = -Inf, xmax = 1994.5, ymin = -Inf, ymax = Inf, alpha = 1, fill="snow2") +
  annotate("rect", xmin = 1994.5, xmax =Inf, ymin = -Inf, ymax = Inf, alpha = 1, fill="snow3")+
  annotate("rect", xmin = 2006.5, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 1, fill="snow4")+
  geom_point(size=3, aes(year,gst,shape=w_trt),inherit.aes = FALSE)+
  scale_shape_manual(name="Rainfall",values=c(1,16)) +
  xlab("Year") +
  theme(plot.margin=unit(c(0,5.5,5.5,5.5),"pt"),legend.position = "bottom",axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())+
  scale_x_continuous(expand = c(0.04,0.04)) +
  coord_cartesian(ylim = c(1,1))+
  guides(color=guide_legend('N deposition',override.aes=list(color=c("snow2","snow3","snow4"),size=5)))+
  geom_vline(xintercept=c(1994.5,2006.5))

lp <- ggplot(subset(simlp, !species %in% "Bromus"),aes(year,median,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  geom_ribbon(aes(x= year, ymin=CIlower,ymax=CIupper, fill=species), colour=NA, alpha = .2) +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"),legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.ticks.x = element_blank())+
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(values=c("#0072B2","#48B99B"),guide=FALSE)+
  scale_fill_manual(name = "Species",labels = c("*Layia*","*Plantago*"),
                    values=c("#0072B2","#48B99B")) +
  scale_x_continuous(expand = c(0.04,0.04))

simulation1 <- ggplot(sim1,aes(year,median,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  geom_ribbon(aes(x= year, ymin=CIlower,ymax=CIupper, fill=species), colour=NA, alpha = .2) +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"),legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_blank()) +
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#48B99B")) +
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                    values=c("#D55E00","#0072B2","#48B99B")) +
  scale_x_continuous(expand = c(0.04,0.04))

simulation2 <- ggplot(sim2,aes(year,median,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  geom_ribbon(aes(x= year, ymin=CIlower,ymax=CIupper, fill=species), color=NA, alpha = .2) +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"),legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_blank()) +
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#48B99B")) +
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                    values=c("#D55E00","#0072B2","#48B99B")) +
  scale_x_continuous(expand = c(0.04,0.04)) #+
 # ylab(expression(atop("Log median",paste("abundance (per  ",m^{2},")")))) 

##legend if needed

leg <- ggplot(sim2,aes(year,median,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  geom_ribbon(aes(x= year, ymin=CIlower,ymax=CIupper, fill=species), alpha = .2) +
  theme(legend.text = element_markdown(),
        strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),
        legend.position = "top",axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#48B99B")) +
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                    values=c("#D55E00","#0072B2","#48B99B")) +
  scale_x_continuous(expand = c(0.04,0.04)) +
  ylab(expression(atop("Log median",paste("abundance (per  ",m^{2},")"))))

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

legend <- as.ggplot(g_legend(leg))

##rainfall dat

ppt <- read.csv(paste(datpath, "JR_rain.csv", sep = "")) %>%
  filter(year %in% 1983:2019)

p2019 <- ggplot(ppt, aes(x=year)) +
  geom_hline(yintercept=565, linetype="dashed",linewidth=0.2)+
  geom_line(aes(y=growing_season_ppt),data=ppt[1:5,],colour="black",size=.8) +
  geom_line(aes(y=growing_season_ppt),data=ppt[5:9,],colour="#bf9b30",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[9:25,],colour="black",size=.8) +
  geom_line(aes(y=growing_season_ppt),data=ppt[25:27,],colour="#bf9b30",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[27:30,],colour="black",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[30:33,],colour="#bf9b30",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[33:37,],colour="black",size=.8)+
  xlab("Year") + 
  #ylab(expression(atop("Rainfall",paste("(mm)"))))+
  ylab("Rainfall (mm)") +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  theme(plot.margin = unit(c(.01,.8,.5,.5), "cm"),axis.text.x = element_text(angle = 90)) +
  annotate("pointrange", x =1983, y =1250.442, 
           ymin = 1250.442, ymax = 1250.442,
           colour = "#593392")+
  annotate("pointrange", x =1998, y =1028.446, 
           ymin = 1028.446, ymax = 1028.446,
           colour = "#593392")+
  annotate("pointrange", x=2017,y=859.028,ymin = 859.028, ymax = 859.028,
           colour = "#593392")+
  scale_y_continuous(limits=c(0,1300),breaks = seq(0,1600,by=565))+
  scale_x_continuous(breaks=seq(1983,2019,by=4))



pdf("sims.pdf",width=6)
plot_grid(lp,simulation2,simulation1,p2019,ncol=1,align="v",rel_heights = c(.015,.015,.015,.015))
dev.off()



