library(tidyverse)
dat <- read.csv("Stan models/parameters.csv")

model.dat <- dat %>%
  pivot_longer(hi.water.hi.N:lo.water.lo.N, names_to ="treatment", values_to="values") %>%
  spread(param, values)

## Set germination and survival fractions from the literature
ps <- .6 # lauren guess
pg <- .92 # gulmon
bs <- .2 # lauren guess
bg <- .98 # gulmon

### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
# new growth function

growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(par.dat$Blambda[i]/(1 + par.dat$BaiP[i]*pg*N$Np[i] + par.dat$BaiB*bg*N$Nb[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(par.dat$Plambda[i]/(1 + par.dat$PaiP[i]*pg*N$Np[i] + par.dat$PaiB[i]*pg*N$Np[i]))
    
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
    myparams <- subset(model.dat, treatment == treatments[trtselect[i]])
    par$Blambda[i] <- as.numeric(subset(myparams, species == "Bromus")[7])
    par$BaiP[i] <- as.numeric(subset(myparams, species == "Bromus")[3])
    par$BaiB[i] <- as.numeric(subset(myparams, species == "Bromus")[6])
    par$Plambda[i] <- as.numeric(subset(myparams, species == "Plantago")[7])
    par$PaiP[i] <- as.numeric(subset(myparams, species == "Plantago")[3])
    par$PaiB[i] <- as.numeric(subset(myparams, species == "Plantago")[6])
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




hi_hi <- grwr(hi_hi, t) %>%
  mutate(treatment = "hi.water.hi.N")
hi_int <- grwr(hi_int, t) %>%
  mutate(treatment = "hi.water.int.N")
hi_lo <- grwr(hi_lo, t) %>%
  mutate(treatment = "hi.water.lo.N")

lo_hi <- grwr(lo_hi, t) %>%
  mutate(treatment = "lo.water.hi.N")
lo_int <- grwr(lo_int, t) %>%
  mutate(treatment = "lo.water.int.N")
lo_lo <- grwr(lo_lo, t) %>%
  mutate(treatment = "lo.water.lo.N")

consistent.out <- rbind(hi_hi, hi_int, hi_lo, lo_hi, lo_int, lo_lo) %>%
  mutate(treatment = as.factor(treatment))


consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment, Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc)) %>%
  separate(treatment, c("water", "delete", "N", "delete2"))

ggplot(consistent.grwr.out, aes(x=N, y=grwrChesson, fill = water)) + geom_bar(stat="identity", position = "dodge") + 
  facet_wrap(~invader)
