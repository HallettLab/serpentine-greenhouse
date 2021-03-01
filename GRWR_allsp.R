library(tidyverse)
library(gridExtra)
dat <- read.csv("Stan models/parameters.csv")

## Data manipulation
model.dat <- dat 
colnames(model.dat)[1] <- "species"
model.dat$species[which(model.dat$species == "Laya")] <- "Layia"
model.dat <- model.dat%>%
  pivot_longer(hi.water.hi.N:lo.water.lo.N, names_to ="treatment", values_to="values") %>%
  spread(param, values) 
colnames(model.dat)[2] <- "treatments"
model.dat[is.na(model.dat)] <- 0

## Set germination and survival fractions from the literature
ps <- .75 # lauren guess
pg <- .92 # gulmon
bs <- .013 # lauren guess
bg <- .98 # gulmon
fs <- .013
fg <- .83
ls <- .15
lg <- .32

###################################
####### Bromus and Plantago #######
###################################

### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
# new growth function

growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(par.dat$Blambda[i]/(1 + par.dat$BaiP[i]*pg*N$Np[i] + par.dat$BaiB*bg*N$Nb[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(par.dat$Plambda[i]/(1 + par.dat$PaiP[i]*pg*N$Np[i] + par.dat$PaiB[i]*bg*N$Nb[i]))
    
  }
  return(N)
}


### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME AND GRWR ####

grwr = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Nb", "Np")
  N1[1,] = c(1,0)
  Bequil <- growth(N1, par.dat, t.num)
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Nb", "Np")
  N2[1,] = c(0,1)
  Pequil <- growth(N2, par.dat, t.num)
  
  N3 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N3) = c("Nb", "Np")
  N3[1,] = c(Bequil[t.num, 1],1)
  Pinvade <- growth(N3, par.dat, t.num)
  Pgrwr <- Pinvade[2,2]/Pinvade[1,2]
  Bgrwc <- Pinvade[2,1]/Pinvade[1,1]
  Pinvade$invader <- "Plantago"
  Pinvade$grwr <- Pgrwr
  Pinvade$time <- as.numeric(row.names(Pinvade))
  Pinvade$Cgrwc <- Bgrwc
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Nb", "Np")
  N4[1,] = c(1, Pequil[t.num, 2])
  Binvade <- growth(N4, par.dat, t.num)
  
  Bgrwr <- Binvade[2,1]/Binvade[1,1]
  Pgrwc <- Binvade[2,2]/Binvade[1,2]
  
  Binvade$invader <- "Bromus"
  Binvade$grwr <- Bgrwr
  Binvade$time <- as.numeric(row.names(Binvade))
  Binvade$Cgrwc <- Pgrwc
  
  
  out <- rbind(Binvade, Pinvade) 
  
  return(out)
}  


### CREATE A FUNCTION THAT SETS CONSISTENT PARAMETERS ACROSS TIMESTEPS ####
consistent.par <- function(model.dat, j.num, t.num){
  par = as.data.frame(matrix(NA, nrow=t, ncol=6))
  colnames(par) = c("Blambda", "BaiP", "BaiB", "Plambda", "PaiP", "PaiB")
  trtselect = rep(j.num, t.num)
  
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatments == treatments[trtselect[i]])
    par$Blambda[i] <- as.numeric(subset(myparams, species == "Bromus")[7])
    par$BaiP[i] <- as.numeric(subset(myparams, species == "Bromus")[6])
    par$BaiB[i] <- as.numeric(subset(myparams, species == "Bromus")[3])
    par$Plambda[i] <- as.numeric(subset(myparams, species == "Plantago")[7])
    par$PaiP[i] <- as.numeric(subset(myparams, species == "Plantago")[6])
    par$PaiB[i] <- as.numeric(subset(myparams, species == "Plantago")[3])
  }
  return(par)
}


### Run the model for the consistent environment ###
t = 100

hi_hi <- consistent.par(model.dat, 1, t)
hi_int <- consistent.par(model.dat, 2, t)
hi_lo <- consistent.par(model.dat, 3, t)
lo_hi <- consistent.par(model.dat, 4, t)
lo_int <- consistent.par(model.dat, 5, t)
lo_lo <- consistent.par(model.dat, 6, t)

hi_hi_grwr <- grwr(hi_hi, t) %>%
  mutate(treatment = "hi.water.hi.N")
hi_int_grwr <- grwr(hi_int, t) %>%
  mutate(treatment = "hi.water.int.N")
hi_lo_grwr <- grwr(hi_lo, t) %>%
  mutate(treatment = "hi.water.lo.N")

lo_hi_grwr <- grwr(lo_hi, t) %>%
  mutate(treatment = "lo.water.hi.N")
lo_int_grwr <- grwr(lo_int, t) %>%
  mutate(treatment = "lo.water.int.N")
lo_lo_grwr <- grwr(lo_lo, t) %>%
  mutate(treatment = "lo.water.lo.N")

consistent.out <- rbind(hi_hi_grwr, hi_int_grwr, hi_lo_grwr, lo_hi_grwr, lo_int_grwr, lo_lo_grwr) %>%
  mutate(treatment = as.factor(treatment))


consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment, Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc)) %>%
  separate(treatment, c("water", "delete", "N", "delete2"))

p1 <- ggplot(consistent.grwr.out, aes(x=N, y=grwrChesson, fill = water)) + geom_bar(stat="identity", position = "dodge") + 
  facet_wrap(~invader)


###################################
####### Bromus and Festuca ########
###################################

### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
# new growth function

growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(par.dat$Blambda[i]/(1 + par.dat$BaiF[i]*fg*N$Nf[i] + par.dat$BaiB*bg*N$Nb[i]))
    
    N$Nf[i+1] = fs*(1-fg)*N$Nf[i] + fg*N$Nf[i]*(par.dat$Flambda[i]/(1 + par.dat$FaiF[i]*fg*N$Nf[i] + par.dat$FaiB[i]*bg*N$Nb[i]))
    
  }
  return(N)
}


### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME AND GRWR ####

grwr = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Nb", "Nf")
  N1[1,] = c(1,0)
  Bequil <- growth(N1, par.dat, t.num)
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Nb", "Nf")
  N2[1,] = c(0,1)
  Fequil <- growth(N2, par.dat, t.num)
  
  N3 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N3) = c("Nb", "Nf")
  N3[1,] = c(Bequil[t.num, 1],1)
  Finvade <- growth(N3, par.dat, t.num)
  Fgrwr <- Finvade[2,2]/Finvade[1,2]
  Bgrwc <- Finvade[2,1]/Finvade[1,1]
  Finvade$invader <- "Festuca"
  Finvade$grwr <- Fgrwr
  Finvade$time <- as.numeric(row.names(Finvade))
  Finvade$Cgrwc <- Bgrwc
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Nb", "Nf")
  N4[1,] = c(1, Fequil[t.num, 2])
  Binvade <- growth(N4, par.dat, t.num)
  
  Bgrwr <- Binvade[2,1]/Binvade[1,1]
  Fgrwc <- Binvade[2,2]/Binvade[1,2]
  
  Binvade$invader <- "Bromus"
  Binvade$grwr <- Bgrwr
  Binvade$time <- as.numeric(row.names(Binvade))
  Binvade$Cgrwc <- Fgrwc
  
  
  out <- rbind(Binvade, Finvade) 
  
  return(out)
}  


### CREATE A FUNCTION THAT SETS CONSISTENT PARAMETERS ACROSS TIMESTEPS ####
consistent.par <- function(model.dat, j.num, t.num){
  par = as.data.frame(matrix(NA, nrow=t, ncol=6))
  colnames(par) = c("Blambda", "BaiP", "BaiB", "Flambda", "FaiF", "FaiB")
  trtselect = rep(j.num, t.num)
  
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatments == treatments[trtselect[i]])
    par$Blambda[i] <- as.numeric(subset(myparams, species == "Bromus")[7])
    par$BaiF[i] <- as.numeric(subset(myparams, species == "Bromus")[4])
    par$BaiB[i] <- as.numeric(subset(myparams, species == "Bromus")[3])
    par$Flambda[i] <- as.numeric(subset(myparams, species == "Festuca")[7])
    par$FaiF[i] <- as.numeric(subset(myparams, species == "Festuca")[4])
    par$FaiB[i] <- as.numeric(subset(myparams, species == "Festuca")[3])
  }
  return(par)
}


### Run the model for the consistent environment ###
t = 100

hi_hi <- consistent.par(model.dat, 1, t)
hi_int <- consistent.par(model.dat, 2, t)
hi_lo <- consistent.par(model.dat, 3, t)
lo_hi <- consistent.par(model.dat, 4, t)
lo_int <- consistent.par(model.dat, 5, t)
lo_lo <- consistent.par(model.dat, 6, t)

hi_hi_grwr <- grwr(hi_hi, t) %>%
  mutate(treatment = "hi.water.hi.N")
hi_int_grwr <- grwr(hi_int, t) %>%
  mutate(treatment = "hi.water.int.N")
hi_lo_grwr <- grwr(hi_lo, t) %>%
  mutate(treatment = "hi.water.lo.N")

lo_hi_grwr <- grwr(lo_hi, t) %>%
  mutate(treatment = "lo.water.hi.N")
lo_int_grwr <- grwr(lo_int, t) %>%
  mutate(treatment = "lo.water.int.N")
lo_lo_grwr <- grwr(lo_lo, t) %>%
  mutate(treatment = "lo.water.lo.N")

consistent.out <- rbind(hi_hi_grwr, hi_int_grwr, hi_lo_grwr, lo_hi_grwr, lo_int_grwr, lo_lo_grwr) %>%
  mutate(treatment = as.factor(treatment))


consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment, Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc)) %>%
  separate(treatment, c("water", "delete", "N", "delete2"))

p2 <- ggplot(consistent.grwr.out, aes(x=N, y=grwrChesson, fill = water)) + geom_bar(stat="identity", position = "dodge") + 
  facet_wrap(~invader)

###################################
####### Bromus and Layia ##########
###################################

### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
# new growth function

growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(par.dat$Blambda[i]/(1 + par.dat$BaiL[i]*lg*N$Nl[i] + par.dat$BaiB*bg*N$Nb[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(par.dat$Llambda[i]/(1 + par.dat$LaiL[i]*lg*N$Nl[i] + par.dat$LaiB[i]*bg*N$Nb[i]))
    
  }
  return(N)
}


### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME AND GRWR ####

grwr = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Nb", "Nl")
  N1[1,] = c(1,0)
  Bequil <- growth(N1, par.dat, t.num)
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Nb", "Nl")
  N2[1,] = c(0,1)
  Lequil <- growth(N2, par.dat, t.num)
  
  N3 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N3) = c("Nb", "Nl")
  N3[1,] = c(Bequil[t.num, 1],1)
  Linvade <- growth(N3, par.dat, t.num)
  Lgrwr <- Linvade[2,2]/Linvade[1,2]
  Bgrwc <- Linvade[2,1]/Linvade[1,1]
  Linvade$invader <- "Layia"
  Linvade$grwr <- Lgrwr
  Linvade$time <- as.numeric(row.names(Linvade))
  Linvade$Cgrwc <- Bgrwc
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Nb", "Nl")
  N4[1,] = c(1, Lequil[t.num, 2])
  Binvade <- growth(N4, par.dat, t.num)
  
  Bgrwr <- Binvade[2,1]/Binvade[1,1]
  Lgrwc <- Binvade[2,2]/Binvade[1,2]
  
  Binvade$invader <- "Bromus"
  Binvade$grwr <- Bgrwr
  Binvade$time <- as.numeric(row.names(Binvade))
  Binvade$Cgrwc <- Lgrwc
  
  
  out <- rbind(Binvade, Linvade) 
  
  return(out)
}  


### CREATE A FUNCTION THAT SETS CONSISTENT PARAMETERS ACROSS TIMESTEPS ####
consistent.par <- function(model.dat, j.num, t.num){
  par = as.data.frame(matrix(NA, nrow=t, ncol=6))
  colnames(par) = c("Blambda", "BaiL", "BaiB", "Llambda", "LaiL", "LaiB")
  trtselect = rep(j.num, t.num)
  
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatments == treatments[trtselect[i]])
    par$Blambda[i] <- as.numeric(subset(myparams, species == "Bromus")[7])
    par$BaiL[i] <- as.numeric(subset(myparams, species == "Bromus")[5])
    par$BaiB[i] <- as.numeric(subset(myparams, species == "Bromus")[3])
    par$Llambda[i] <- as.numeric(subset(myparams, species == "Layia")[7])
    par$LaiL[i] <- as.numeric(subset(myparams, species == "Layia")[5])
    par$LaiB[i] <- as.numeric(subset(myparams, species == "Layia")[3])
  }
  return(par)
}


### Run the model for the consistent environment ###
t = 100

hi_hi <- consistent.par(model.dat, 1, t)
hi_int <- consistent.par(model.dat, 2, t)
hi_lo <- consistent.par(model.dat, 3, t)
lo_hi <- consistent.par(model.dat, 4, t)
lo_int <- consistent.par(model.dat, 5, t)
lo_lo <- consistent.par(model.dat, 6, t)

hi_hi_grwr <- grwr(hi_hi, t) %>%
  mutate(treatment = "hi.water.hi.N")
hi_int_grwr <- grwr(hi_int, t) %>%
  mutate(treatment = "hi.water.int.N")
hi_lo_grwr <- grwr(hi_lo, t) %>%
  mutate(treatment = "hi.water.lo.N")

lo_hi_grwr <- grwr(lo_hi, t) %>%
  mutate(treatment = "lo.water.hi.N")
lo_int_grwr <- grwr(lo_int, t) %>%
  mutate(treatment = "lo.water.int.N")
lo_lo_grwr <- grwr(lo_lo, t) %>%
  mutate(treatment = "lo.water.lo.N")

consistent.out <- rbind(hi_hi_grwr, hi_int_grwr, hi_lo_grwr, lo_hi_grwr, lo_int_grwr, lo_lo_grwr) %>%
  mutate(treatment = as.factor(treatment))


consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment, Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc)) %>%
  separate(treatment, c("water", "delete", "N", "delete2"))

p3 <- ggplot(consistent.grwr.out, aes(x=N, y=grwrChesson, fill = water)) + geom_bar(stat="identity", position = "dodge") + 
  facet_wrap(~invader)


###################################
####### Festuca and Plantago ######
###################################

### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
# new growth function

growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Nf[i+1] = fs*(1-fg)*N$Nf[i] + fg*N$Nf[i]*(par.dat$Flambda[i]/(1 + par.dat$FaiF[i]*fg*N$Nf[i] + par.dat$FaiP[i]*pg*N$Np[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(par.dat$Plambda[i]/(1 + par.dat$PaiP[i]*pg*N$Np[i] + par.dat$PaiF[i]*fg*N$Nf[i]))
    
  }
  return(N)
}


### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME AND GRWR ####

grwr = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Nf", "Np")
  N1[1,] = c(1,0)
  Fequil <- growth(N1, par.dat, t.num)
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Nf", "Np")
  N2[1,] = c(0,1)
  Pequil <- growth(N2, par.dat, t.num)
  
  N3 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N3) = c("Nf", "Np")
  N3[1,] = c(Fequil[t.num, 1],1)
  Pinvade <- growth(N3, par.dat, t.num)
  Pgrwr <- Pinvade[2,2]/Pinvade[1,2]
  Fgrwc <- Pinvade[2,1]/Pinvade[1,1]
  Pinvade$invader <- "Plantago"
  Pinvade$grwr <- Pgrwr
  Pinvade$time <- as.numeric(row.names(Pinvade))
  Pinvade$Cgrwc <- Fgrwc
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Nf", "Np")
  N4[1,] = c(1, Pequil[t.num, 2])
  Finvade <- growth(N4, par.dat, t.num)
  
  Fgrwr <- Finvade[2,1]/Finvade[1,1]
  Pgrwc <- Finvade[2,2]/Finvade[1,2]
  
  Finvade$invader <- "Festuca"
  Finvade$grwr <- Fgrwr
  Finvade$time <- as.numeric(row.names(Finvade))
  Finvade$Cgrwc <- Pgrwc
  
  
  out <- rbind(Finvade, Pinvade) 
  
  return(out)
}  


### CREATE A FUNCTION THAT SETS CONSISTENT PARAMETERS ACROSS TIMESTEPS ####
consistent.par <- function(model.dat, j.num, t.num){
  par = as.data.frame(matrix(NA, nrow=t, ncol=6))
  colnames(par) = c("Flambda", "FaiP", "FaiF", "Plambda", "PaiP", "PaiF")
  trtselect = rep(j.num, t.num)
  
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatments == treatments[trtselect[i]])
    par$Flambda[i] <- as.numeric(subset(myparams, species == "Festuca")[7])
    par$FaiP[i] <- as.numeric(subset(myparams, species == "Festuca")[6])
    par$FaiF[i] <- as.numeric(subset(myparams, species == "Festuca")[4])
    par$Plambda[i] <- as.numeric(subset(myparams, species == "Plantago")[7])
    par$PaiP[i] <- as.numeric(subset(myparams, species == "Plantago")[6])
    par$PaiF[i] <- as.numeric(subset(myparams, species == "Plantago")[4])
  }
  return(par)
}


### Run the model for the consistent environment ###
t = 100

hi_hi <- consistent.par(model.dat, 1, t)
hi_int <- consistent.par(model.dat, 2, t)
hi_lo <- consistent.par(model.dat, 3, t)
lo_hi <- consistent.par(model.dat, 4, t)
lo_int <- consistent.par(model.dat, 5, t)
lo_lo <- consistent.par(model.dat, 6, t)

hi_hi_grwr <- grwr(hi_hi, t) %>%
  mutate(treatment = "hi.water.hi.N")
hi_int_grwr <- grwr(hi_int, t) %>%
  mutate(treatment = "hi.water.int.N")
hi_lo_grwr <- grwr(hi_lo, t) %>%
  mutate(treatment = "hi.water.lo.N")

lo_hi_grwr <- grwr(lo_hi, t) %>%
  mutate(treatment = "lo.water.hi.N")
lo_int_grwr <- grwr(lo_int, t) %>%
  mutate(treatment = "lo.water.int.N")
lo_lo_grwr <- grwr(lo_lo, t) %>%
  mutate(treatment = "lo.water.lo.N")

consistent.out <- rbind(hi_hi_grwr, hi_int_grwr, hi_lo_grwr, lo_hi_grwr, lo_int_grwr, lo_lo_grwr) %>%
  mutate(treatment = as.factor(treatment))


consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment, Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc)) %>%
  separate(treatment, c("water", "delete", "N", "delete2"))

p4 <- ggplot(consistent.grwr.out, aes(x=N, y=grwrChesson, fill = water)) + geom_bar(stat="identity", position = "dodge") + 
  facet_wrap(~invader)


###################################
####### Layia and Plantago ######
###################################

### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
# new growth function

growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(par.dat$Llambda[i]/(1 + par.dat$LaiL[i]*lg*N$Nl[i] + par.dat$LaiP[i]*pg*N$Np[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(par.dat$Plambda[i]/(1 + par.dat$PaiP[i]*pg*N$Np[i] + par.dat$PaiL[i]*lg*N$Nl[i]))
    
  }
  return(N)
}


### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME AND GRWR ####

grwr = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Nl", "Np")
  N1[1,] = c(1,0)
  Lequil <- growth(N1, par.dat, t.num)
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Nl", "Np")
  N2[1,] = c(0,1)
  Pequil <- growth(N2, par.dat, t.num)
  
  N3 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N3) = c("Nl", "Np")
  N3[1,] = c(Lequil[t.num, 1],1)
  Pinvade <- growth(N3, par.dat, t.num)
  Pgrwr <- Pinvade[2,2]/Pinvade[1,2]
  Lgrwc <- Pinvade[2,1]/Pinvade[1,1]
  Pinvade$invader <- "Plantago"
  Pinvade$grwr <- Pgrwr
  Pinvade$time <- as.numeric(row.names(Pinvade))
  Pinvade$Cgrwc <- Lgrwc
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Nl", "Np")
  N4[1,] = c(1, Pequil[t.num, 2])
  Linvade <- growth(N4, par.dat, t.num)
  
  Lgrwr <- Linvade[2,1]/Linvade[1,1]
  Pgrwc <- Linvade[2,2]/Linvade[1,2]
  
  Linvade$invader <- "Layia"
  Linvade$grwr <- Lgrwr
  Linvade$time <- as.numeric(row.names(Linvade))
  Linvade$Cgrwc <- Pgrwc
  
  
  out <- rbind(Linvade, Pinvade) 
  
  return(out)
}  


### CREATE A FUNCTION THAT SETS CONSISTENT PARAMETERS ACROSS TIMESTEPS ####
consistent.par <- function(model.dat, j.num, t.num){
  par = as.data.frame(matrix(NA, nrow=t, ncol=6))
  colnames(par) = c("Llambda", "LaiP", "LaiL", "Plambda", "PaiP", "PaiL")
  trtselect = rep(j.num, t.num)
  
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatments == treatments[trtselect[i]])
    par$Llambda[i] <- as.numeric(subset(myparams, species == "Layia")[7])
    par$LaiP[i] <- as.numeric(subset(myparams, species == "Layia")[6])
    par$LaiL[i] <- as.numeric(subset(myparams, species == "Layia")[5])
    par$Plambda[i] <- as.numeric(subset(myparams, species == "Plantago")[7])
    par$PaiP[i] <- as.numeric(subset(myparams, species == "Plantago")[6])
    par$PaiL[i] <- as.numeric(subset(myparams, species == "Plantago")[5])
  }
  return(par)
}


### Run the model for the consistent environment ###
t = 100

hi_hi <- consistent.par(model.dat, 1, t)
hi_int <- consistent.par(model.dat, 2, t)
hi_lo <- consistent.par(model.dat, 3, t)
lo_hi <- consistent.par(model.dat, 4, t)
lo_int <- consistent.par(model.dat, 5, t)
lo_lo <- consistent.par(model.dat, 6, t)

hi_hi_grwr <- grwr(hi_hi, t) %>%
  mutate(treatment = "hi.water.hi.N")
hi_int_grwr <- grwr(hi_int, t) %>%
  mutate(treatment = "hi.water.int.N")
hi_lo_grwr <- grwr(hi_lo, t) %>%
  mutate(treatment = "hi.water.lo.N")

lo_hi_grwr <- grwr(lo_hi, t) %>%
  mutate(treatment = "lo.water.hi.N")
lo_int_grwr <- grwr(lo_int, t) %>%
  mutate(treatment = "lo.water.int.N")
lo_lo_grwr <- grwr(lo_lo, t) %>%
  mutate(treatment = "lo.water.lo.N")

consistent.out <- rbind(hi_hi_grwr, hi_int_grwr, hi_lo_grwr, lo_hi_grwr, lo_int_grwr, lo_lo_grwr) %>%
  mutate(treatment = as.factor(treatment))


consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment, Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc)) %>%
  separate(treatment, c("water", "delete", "N", "delete2"))

p5 <- ggplot(consistent.grwr.out, aes(x=N, y=grwrChesson, fill = water)) + geom_bar(stat="identity", position = "dodge") + 
  facet_wrap(~invader)


###################################
####### Layia and Festuca ######
###################################

### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
# new growth function

growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(par.dat$Llambda[i]/(1 + par.dat$LaiL[i]*lg*N$Nl[i] + par.dat$LaiF[i]*fg*N$Nf[i]))
    
    N$Nf[i+1] = fs*(1-fg)*N$Nf[i] + fg*N$Nf[i]*(par.dat$Flambda[i]/(1 + par.dat$FaiF[i]*fg*N$Nf[i] + par.dat$FaiL[i]*lg*N$Nl[i]))
    
  }
  return(N)
}


### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME AND GRWR ####

grwr = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Nl", "Nf")
  N1[1,] = c(1,0)
  Lequil <- growth(N1, par.dat, t.num)
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Nl", "Nf")
  N2[1,] = c(0,1)
  Fequil <- growth(N2, par.dat, t.num)
  
  N3 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N3) = c("Nl", "Nf")
  N3[1,] = c(Lequil[t.num, 1],1)
  Finvade <- growth(N3, par.dat, t.num)
  Fgrwr <- Finvade[2,2]/Finvade[1,2]
  Lgrwc <- Finvade[2,1]/Finvade[1,1]
  Finvade$invader <- "Festuca"
  Finvade$grwr <- Fgrwr
  Finvade$time <- as.numeric(row.names(Finvade))
  Finvade$Cgrwc <- Lgrwc
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Nl", "Nf")
  N4[1,] = c(1, Fequil[t.num, 2])
  Linvade <- growth(N4, par.dat, t.num)
  
  Lgrwr <- Linvade[2,1]/Linvade[1,1]
  Fgrwc <- Linvade[2,2]/Linvade[1,2]
  
  Linvade$invader <- "Layia"
  Linvade$grwr <- Lgrwr
  Linvade$time <- as.numeric(row.names(Linvade))
  Linvade$Cgrwc <- Fgrwc
  
  
  out <- rbind(Linvade, Finvade) 
  
  return(out)
}  


### CREATE A FUNCTION THAT SETS CONSISTENT PARAMETERS ACROSS TIMESTEPS ####
consistent.par <- function(model.dat, j.num, t.num){
  par = as.data.frame(matrix(NA, nrow=t, ncol=6))
  colnames(par) = c("Llambda", "LaiP", "LaiL", "Flambda", "FaiF", "FaiL")
  trtselect = rep(j.num, t.num)
  
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatments == treatments[trtselect[i]])
    par$Llambda[i] <- as.numeric(subset(myparams, species == "Layia")[7])
    par$LaiF[i] <- as.numeric(subset(myparams, species == "Layia")[4])
    par$LaiL[i] <- as.numeric(subset(myparams, species == "Layia")[5])
    par$Flambda[i] <- as.numeric(subset(myparams, species == "Festuca")[7])
    par$FaiF[i] <- as.numeric(subset(myparams, species == "Festuca")[4])
    par$FaiL[i] <- as.numeric(subset(myparams, species == "Festuca")[5])
  }
  return(par)
}


### Run the model for the consistent environment ###
t = 100

hi_hi <- consistent.par(model.dat, 1, t)
hi_int <- consistent.par(model.dat, 2, t)
hi_lo <- consistent.par(model.dat, 3, t)
lo_hi <- consistent.par(model.dat, 4, t)
lo_int <- consistent.par(model.dat, 5, t)
lo_lo <- consistent.par(model.dat, 6, t)

hi_hi_grwr <- grwr(hi_hi, t) %>%
  mutate(treatment = "hi.water.hi.N")
hi_int_grwr <- grwr(hi_int, t) %>%
  mutate(treatment = "hi.water.int.N")
hi_lo_grwr <- grwr(hi_lo, t) %>%
  mutate(treatment = "hi.water.lo.N")

lo_hi_grwr <- grwr(lo_hi, t) %>%
  mutate(treatment = "lo.water.hi.N")
lo_int_grwr <- grwr(lo_int, t) %>%
  mutate(treatment = "lo.water.int.N")
lo_lo_grwr <- grwr(lo_lo, t) %>%
  mutate(treatment = "lo.water.lo.N")

consistent.out <- rbind(hi_hi_grwr, hi_int_grwr, hi_lo_grwr, lo_hi_grwr, lo_int_grwr, lo_lo_grwr) %>%
  mutate(treatment = as.factor(treatment))


consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment, Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc)) %>%
  separate(treatment, c("water", "delete", "N", "delete2"))

p6 <- ggplot(consistent.grwr.out, aes(x=N, y=grwrChesson, fill = water)) + geom_bar(stat="identity", position = "dodge") + 
  facet_wrap(~invader)
