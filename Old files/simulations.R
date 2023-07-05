library(tidyverse)
library(ggtext)
library(cowplot)
library(grid)
library(ggplotify)

## Read in data
params <- read.csv(paste(datpath, "params.csv", sep = "")) # parameters from first stan model fits
trt <- read.csv(paste(datpath, "years_trt.csv", sep = "")) %>%
  rename(w_trt=type_year)
cover <-read.csv(paste(datpath,"JR_cover_1mplot.csv",sep=""))

## Determine equilibrium conditions for each species in isolation
params_dat <- params %>%
  separate(treatments, into = c("w_trt","n_trt"),sep=-5)%>%
  filter(species == "Bromus" |species=="Plantago"|species=="Layia")
params_dat$w_trt[params_dat$w_trt == c("hi.water","hi.water.","hi.water")] <- "Wet"
params_dat$w_trt[params_dat$w_trt == c("lo.water", "lo.water.","lo.water")] <- "Dry"
params_dat$n_trt[params_dat$n_trt == ".hi.N"] <- "hi.N"
params_dat$n_trt[params_dat$n_trt == ".lo.N"] <- "lo.N"
params_dat$w_trt <- as.factor(params_dat$w_trt)
params_dat$n_trt <- as.factor(params_dat$n_trt)

pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0*g)
  return(N)
}

## Bromus
bromus <- params_dat %>%
  filter(species == "Bromus")

N0 <- 27 #average number of seeds by weight added to greenhouse pots
time <- length(trt$w_trt)
N_bromus <- rep(NA, time)
N_bromus[1] <- N0
bs <- 0.013
bg <- 0.98

for (t in 1:time) {
  params <- subset(bromus,w_trt==trt$w_trt[t] & n_trt==trt$n_trt[t])
  N_bromus[t+1] <- pop.equilibrium(N0=N_bromus[t], s=bs, g=bg, a_intra=bromus$alpha_bromus, lambda=bromus$lambda)
}

# check output
plot(seq(1:(time+1)), N_bromus, type="l")
N_bromus

## Layia
layia <- params_dat %>%
  filter(species == "Layia")

N0 <- 89 #average number of seeds by weight added to greenhouse pots
time <- length(trt$w_trt)
N_layia <- rep(NA, time)
N_layia[1] <- N0
ls <- 0.15
lg <- 0.32

for (t in 1:time) {
  params <- subset(layia,w_trt==trt$w_trt[t] & n_trt==trt$n_trt[t])
  N_layia[t+1] <- pop.equilibrium(N0=N_layia[t], s=ls, g=lg, a_intra=layia$alpha_lapl, lambda=layia$lambda)
}

# check output
plot(seq(1:(time+1)), N_layia, type="l")
N_layia

## Plantago
pler <- params_dat %>%
  filter(species == "Plantago")

N0 <- 26 #average number of seeds by weight added to greenhouse pots
time <- length(trt$w_trt)
N_pler <- rep(NA, time)
N_pler[1] <- N0
ps <- 0.92
pg <- 0.75

for (t in 1:time) {
  params <- subset(pler,w_trt==trt$w_trt[t] & n_trt==trt$n_trt[t])
  N_pler[t+1] <- pop.equilibrium(N0=N_pler[t], s=ps, g=pg, a_intra=pler$alpha_plantago, lambda=pler$lambda)
}

# check output
plot(seq(1:(time+1)), N_pler, type="l")
N_pler


## Data manipulation for starting conditions
# Equilibrium abundances
brho_eq <- N_bromus[37]
pler_eq <- N_pler[37]
lapl_eq <- N_layia[37]

cover_dat <- cover %>%
  filter(species == "BRMO" | species == "PLER" | species == "LAPL") %>%
  group_by(year, species) %>%
  summarize(mean_cover = mean(cover)) %>%
  filter(year == 1983) %>%
  mutate(species = c("Bromus","Layia","Plantago")) %>%
  mutate(abundance = mean_cover*c(brho_eq,lapl_eq,pler_eq)/100)%>%
  select(-mean_cover) %>%
  pivot_wider(names_from = species, values_from = abundance)

params2_dat <- read.csv(paste(datpath, "params2.csv", sep = ""))%>%
  separate(treatments, into = c("w_trt","n_trt"),sep=-5)%>%
  select(-4,-9,-14,-18,-19,-20,-21,-22)
params2_dat$w_trt[params2_dat$w_trt == c("hi.water","hi.water.","hi.water")] <- "Wet"
params2_dat$w_trt[params2_dat$w_trt == c("lo.water", "lo.water.","lo.water")] <- "Dry"
params2_dat$n_trt[params2_dat$n_trt == ".hi.N"] <- "hi.N"
params2_dat$n_trt[params2_dat$n_trt == ".lo.N"] <- "lo.N"
params2_dat$w_trt <- as.factor(params2_dat$w_trt)
params2_dat$n_trt <- as.factor(params2_dat$n_trt) # parameters in a different format

cov_trt <- left_join(cover_dat,trt)
trt_params <- left_join(trt,params2_dat)

start_dat <- full_join(cov_trt, trt_params) 

#####################################
########## Simulation ###############
#####################################

start_dat <- start_dat %>%
  filter(year !=1983) 

year <- length(1983:2019)

## all species sim 1 without altering bromus germination
N = as.data.frame(matrix(NA, nrow=37, ncol=3))
colnames(N) = c("Nb", "Nl","Np")
N[1,] = c(cov_trt$Bromus,cov_trt$Layia,cov_trt$Plantago)

growth = function(N, start_dat,year){
  for (i in 1:(year-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(start_dat$blambda[i]/(1 + start_dat$bap[i]*pg*N$Np[i] + start_dat$bab[i]*bg*N$Nb[i] + start_dat$bal[i]*lg*N$Nl[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pab[i]*bg*N$Nb[i] + start_dat$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$lab[i]*bg*N$Nb[i] + start_dat$lap[i]*pg*N$Np[i]))
   
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
    
  }

  return(N)
}

Nblp <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")
Nblp$species[Nblp$species == "Nb"] <- "Bromus"
Nblp$species[Nblp$species == "Nl"] <- "Layia"
Nblp$species[Nblp$species == "Np"] <- "Plantago"

blp <- ggplot(Nblp,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3)  +
  theme(legend.text = element_markdown(),legend.position = "top",axis.title.x = element_blank(),axis.text.x = element_blank()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#009E73"))

## all species sim with bromus germination = 0
#BRHO, PLER, and LAPL
brho_bg <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98)
start_dat_brho <- start_dat %>%
  cbind(bg = brho_bg)

N = as.data.frame(matrix(NA, nrow=37, ncol=3))
colnames(N) = c("Nb", "Np","Nl")
N[1,] =c(cov_trt$Bromus,cov_trt$Layia,cov_trt$Plantago)

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

Nblp2 <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")
Nblp2$species[Nblp2$species == "Nb"] <- "Bromus"
Nblp2$species[Nblp2$species == "Np"] <- "Plantago"
Nblp2$species[Nblp2$species == "Nl"] <- "Layia"

blp2 <- ggplot(Nblp2,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3) + xlab("Year") +
  theme(legend.position= "none") +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                     values=c("#D55E00","#0072B2","#009E73"))

## PLER and LAPL sim
N = as.data.frame(matrix(NA, nrow=37, ncol=2))
colnames(N) = c("Np", "Nl")
N[1,] =c(cov_trt$Plantago,cov_trt$Layia)


growth = function(N, start_dat,year){
  for (i in 1:(year-1)){
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$lap[i]*pg*N$Np[i]))
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
  }
  
  return(N)
}

Nlp <- growth(N,start_dat,year) %>%
  mutate(year = 1983:2019) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")
Nlp$species[Nlp$species == "Np"] <- "Plantago"
Nlp$species[Nlp$species == "Nl"] <- "Layia"

lp <- ggplot(Nlp,aes(year,abundance,color=species)) + geom_line(size = .8) + geom_point(size=1.3)  +
  theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank()) +
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(name = "Species",labels = c("*Layia*","*Plantago*"),
                     values=c("#0072B2","#009E73"))


#########################
#####Visualization#######
#########################
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 16),
              axis.text = element_text(size = 16))

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
  geom_point(size=3, aes(year,gst,shape=w_trt),inherit.aes = FALSE)+
  scale_shape_manual(name="Rainfall",values=c(1,16)) +
  xlab("Year") +
  theme(plot.margin=unit(c(0,5.5,5.5,5.5),"pt"),legend.position = "bottom",axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())+
  scale_x_continuous(expand = c(0.04,0.04)) +
  coord_cartesian(ylim = c(1,1))+
  guides(color=guide_legend('N deposition',override.aes=list(color=c("snow2","snow3","snow4"),size=5)))+
  geom_vline(xintercept=c(1994.5,2006.5))

plot_grid(legend,lp,sim1,sim2,bar,ncol=1,align="v",rel_heights = c(.2,1,1,1,.6),labels = c("","a)","b)","c)","d)"))


