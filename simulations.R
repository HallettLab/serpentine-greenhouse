library(tidyverse)
library(ggtext)
library(cowplot)
library(grid)
library(ggplotify)

## Read in data
params2 <- read.csv(paste(datpath, "params2.csv", sep = ""))
trt <- read.csv(paste(datpath, "years_trt.csv", sep = "")) 
cover <-read.csv(paste(datpath,"JR_cover_1mplot.csv",sep=""))

## Determine equilibrium conditions for each species in isolation
params_dat <- params2

pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0*g)
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

lp <- ggplot(N,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3)  +
  theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Layia*","*Plantago*"),
                     values=c("#0072B2","#009E73"))
Nlp <- N

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

blp <- ggplot(N,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3)  +
  theme(legend.text = element_markdown(),legend.position = "top",axis.title.x = element_blank(),axis.text.x = element_blank()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#009E73"))
Nblp <- N

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
  theme(legend.position= "none") +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#009E73"))
Nblp2 <- N

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

ggarrange(blp,blp2,
          labels = c("a)", "b)"),
          ncol = 1, nrow = 2)

#two plots of simulations
Nblp <- Nblp %>%
  inner_join(trt)

Nblp2 <- Nblp2 %>%
  inner_join(trt)

Nlp <- Nlp %>%
  inner_join(trt)


blpN <- inner_join(Nblp,Nblp2) 

blp2N <- blpN%>%
  pivot_longer(3:4,names_to = "sim",values_to = "abundance") %>%
  add_column(N = if_else(.$year < 1995,"Low",ifelse(.$year > 1995 & .$year > 2006, "High", "Intermediate")))

blp2N$N <- factor(blp2N$N, levels = c("Low","Intermediate","High"))



p <- ggplot(blp2N,aes(year,abundance,color=species)) +
  annotate("rect", xmin = 1983, xmax = 1995, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow2")+
  annotate("rect", xmin = 1995, xmax = 2007, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow3")+
  annotate("rect", xmin = 2007, xmax = 2019, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow4")+
  geom_line(size = .8) + geom_point(size=1.3) + xlab("Year") +
  theme(legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "top") +
  facet_wrap(~sim,ncol=1)+
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#009E73")) 

trts <- ggplot(trt,aes(year,growing_season_ppt)) +
  annotate("rect", xmin = -Inf, xmax = 1995, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow2") +
  annotate("rect", xmin = 1995, xmax = 2007, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow3")+
  annotate("rect", xmin = 2007, xmax = Inf, min = -0.5, ymax = Inf, alpha = 0.6, fill="snow4")+
  geom_line(size=0.8)+
  geom_hline(yintercept=565,lty=2)+
  geom_point(size=3, aes(shape=type_year))+
  scale_shape_manual(name="Growing season type",values=c(1,16)) +
  xlab("Year") +
  ylab("Precipitation (mm)")+
  theme(legend.position = "top",axis.text.x = element_blank(),axis.title.x = element_blank(),
        plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"))+
  scale_x_continuous(expand = c(0.04,0.04)) +
  annotate("text", x=1988,y=1200, label="Low N")+
  annotate("text", x=2001,y=1200, label="Intermediate N") +
  annotate("text", x=2014,y=1200, label="High N") 
  


lp <- ggplot(Nlp,aes(year,abundance,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"),legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(values=c("#0072B2","#009E73"),guide=FALSE)+
  scale_x_continuous(expand = c(0.04,0.04))



sim1 <- ggplot(Nblp,aes(year,abundance,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"),legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#009E73")) +
  scale_x_continuous(expand = c(0.04,0.04))

leg <- ggplot(Nblp,aes(year,abundance,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "top",axis.text.x = element_blank(),axis.title.x = element_blank()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#009E73")) +
  scale_x_continuous(expand = c(0.04,0.04))+
  theme(legend.text = element_markdown())

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

legend <- as.ggplot(g_legend(leg))



sim2 <-  ggplot(Nblp2,aes(year,abundance,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  theme(plot.margin=unit(c(5.5,10,0,5.5),"pt"),legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "none",axis.title.y = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#009E73")) +
  scale_x_continuous(expand = c(0.04, 0.04))


trt2 <- trt %>%
  mutate(gst = 1,n=NA) 

trt2$n[trt2$n_trt == "lo.N"] <- "Low"
trt2$n[trt2$n_trt == "int.N"] <- "Intermediate"
trt2$n[trt2$n_trt == "hi.N"] <- "High"

trt2$n <- factor(trt2$n,levels = c("Low","Intermediate","High"))

bar <- ggplot(trt2,aes(year,gst)) +
  geom_point(aes(year,gst,color=n),size=0,shape =15)+
  annotate("rect", xmin = -Inf, xmax = 1994.5, ymin = -Inf, ymax = Inf, alpha = 1, fill="snow2") +
  annotate("rect", xmin = 1994.5, xmax = 2006.5, ymin = -Inf, ymax = Inf, alpha = 1, fill="snow3")+
  annotate("rect", xmin = 2006.5, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 1, fill="snow4")+
  geom_point(size=3, aes(year,gst,shape=type_year),inherit.aes = FALSE)+
  scale_shape_manual(name="Rainfall",values=c(1,16)) +
  xlab("Year") +
  theme(plot.margin=unit(c(0,5.5,5.5,5.5),"pt"),legend.position = "bottom",axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())+
  scale_x_continuous(expand = c(0.04,0.04)) +
  coord_cartesian(ylim = c(1,1))+
  guides(color=guide_legend('N deposition',override.aes=list(color=c("snow2","snow3","snow4"),size=5)))+
  geom_vline(xintercept=c(1994.5,2006.5))


plot_grid(legend,lp,sim1,sim2,bar,ncol=1,align="v",rel_heights = c(.2,1,1,1,.6),labels = c("","a)","b)","c)","d)"))


