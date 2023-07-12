library(tidyverse)
library(ggtext)
library(cowplot)
library(grid)
library(ggplotify)
library(rstan)
library(here)

ps <- .75 # gulmon
pg <- .92 # gulmon
bs <- .013 # andrew
bg <- .98 # gulmon
ls <- .15 # rossington
lg <- .32 # rossington

## Read in data
#params2 <- read.csv(paste(datpath, "params2_LGS.csv", sep = "")) # parameters from first stan model fits
trt <- read.csv(paste(datpath, "years_trt.csv", sep = "")) %>%
  rename(w_trt=type_year)
cover <-read.csv(paste(datpath,"JR_cover_1mplot.csv",sep=""))
stems_dat <-  read.csv(paste(datpath, "/stems_background.csv", sep = "")) %>%
  filter(seed_sp != "FEMI") %>%
  mutate(recruit=stem_density/seed_added) 

##survival and germination fractions
ps <- .75 # gulmon
pg <- .92 # gulmon
bs <- .013 # andrew
bg <- .98 # gulmon
ls <- .15 # rossington
lg <- .32 # rossington

##Load stan models from Stan models folder
load(here("Stan models","no_dist_seeds_pler_hi_hi.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_hi.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_hi_hi.Rdata"))
load(here("Stan models","no_dist_seeds_pler_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_hi_int.Rdata"))
load(here("Stan models","no_dist_seeds_pler_hi_lo.Rdata"))
load(here("Stan models","no_dist_seeds_brho_hi_lo.Rdata"))
load(here("Stan models","no_dist_seeds_pler_lo_hi.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_hi.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_lo_hi.Rdata"))
load(here("Stan models","no_dist_seeds_pler_lo_int.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_int.Rdata"))
load(here("Stan models","no_dist_seeds_lapl_lo_int.Rdata"))
load(here("Stan models","no_dist_seeds_pler_lo_lo.Rdata"))
load(here("Stan models","no_dist_seeds_brho_lo_lo.Rdata"))

#Extract parameters from models and rename 
brho_hi_hi <- rstan::extract(no_dist_seeds_brho_hi_hi)
brho_hi_int <- rstan::extract(no_dist_seeds_brho_hi_int)
brho_hi_lo <- rstan::extract(no_dist_seeds_brho_hi_lo)
brho_lo_hi <- rstan::extract(no_dist_seeds_brho_lo_hi)
brho_lo_int <- rstan::extract(no_dist_seeds_brho_lo_int)
brho_lo_lo <- rstan::extract(no_dist_seeds_brho_lo_lo)
lapl_hi_hi <- rstan::extract(no_dist_seeds_lapl_hi_hi)
lapl_hi_int <- rstan::extract(no_dist_seeds_lapl_hi_int)
lapl_lo_hi <- rstan::extract(no_dist_seeds_lapl_lo_hi)
lapl_lo_int <- rstan::extract(no_dist_seeds_lapl_lo_int)
pler_hi_hi <- rstan::extract(no_dist_seeds_pler_hi_hi)
pler_hi_int <- rstan::extract(no_dist_seeds_pler_hi_int)
pler_hi_lo <- rstan::extract(no_dist_seeds_pler_hi_lo)
pler_lo_hi <- rstan::extract(no_dist_seeds_pler_lo_hi)
pler_lo_int <- rstan::extract(no_dist_seeds_pler_lo_int)
pler_lo_lo <- rstan::extract(no_dist_seeds_pler_lo_lo)

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

##join posts with trt dataframe
cover_dat <- cover %>%
  filter(treatment == "c") %>%
  filter(species == "PLER" | species == "BRMO" | species == "LAPL") %>%
  filter(year==1983) %>%
  group_by(species) %>%
  summarize(mean_cover=mean(cover)/100) %>%
  pivot_wider(names_from=species,values_from=mean_cover) %>%
  mutate(year=1983)

cover_trt <- left_join(cover_dat,trt) %>%
  rename(BRHO = BRMO) %>%
  mutate(BRHO = BRHO*16,PLER=PLER*128,LAPL=LAPL*0)

start_dat <- left_join(trt_posts,cover_trt) 

#####################################
########## Simulation ###############
#####################################

start_dat <- start_dat %>%
  filter(year !=1983) 

#year <- length(1:72001)

## all species sim 1 without altering bromus germination
N = as.data.frame(matrix(NA, nrow=72001, ncol=3))
colnames(N) = c("Nb", "Nl","Np")
N[1,] = c(cover_trt$BRHO,cover_trt$LAPL,cover_trt$PLER)

growth = function(N, start_dat,replicate.var){
  replicate.var = N[1,1]
  for (i in 1:(replicate.var-1)){
    N$Nb[i+1] = bs*(1-bg)*N$Nb[i]  + bg*N$Nb[i]*(start_dat$blambda[i]/(1 + start_dat$bap[i]*pg*N$Np[i] + start_dat$bab[i]*bg*N$Nb[i] + start_dat$bal[i]*lg*N$Nl[i]))
    
    N$Np[i+1] = ps*(1-pg)*N$Np[i] + pg*N$Np[i]*(start_dat$plambda[i]/(1 + start_dat$pap[i]*pg*N$Np[i] + start_dat$pab[i]*bg*N$Nb[i] + start_dat$pal[i]*lg*N$Nl[i]))
    
    N$Nl[i+1] = ls*(1-lg)*N$Nl[i] + lg*N$Nl[i]*(start_dat$llambda[i]/(1 + start_dat$lal[i]*lg*N$Nl[i] + start_dat$lab[i]*bg*N$Nb[i] + start_dat$lap[i]*pg*N$Np[i]))
   
    N$Nl[i+1] = ifelse(N$Nl[i+1]<0.1,0.1,N$Nl[i+1])
    
    
  }

  return(N)
}

X[[1]]


X[[1]]


index <- cbind(replicate.var=rep(seq(1:2000),36))
start <- cbind(start_dat,index) 
start <- start[order(start[["replicate.var"]]),]
X <- split(start, start["replicate.var"])
out <- lapply(X, FUN=growth, start_dat)
reps <- unique(start["replicate.var"])
output <- cbind(reps, do.call("rbind", out))
names(output) = c("replicate.var", "year","water_trt","n_trt","blambda","bab",
                  "bap","bal","llambda","lal","lap","lab","plambda","pap",
                  "pab","pal","BRHO","LAPL","PLER")

growth(X[1], replicate.var = "year")

Nblp <- growth(N,start_dat,year) %>%
  pivot_longer(!year,names_to="species",values_to ="abundance")%>%
Nblp$species[Nblp$species == "Nb"] <- "Bromus"
Nblp$species[Nblp$species == "Nl"] <- "Layia"
Nblp$species[Nblp$species == "Np"] <- "Plantago"

## all species sim with bromus germination = 0
#BRHO, PLER, and LAPL
brho_bg <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98)
start_dat_brho <- start_dat %>%
  cbind(bg = brho_bg)

N = as.data.frame(matrix(NA, nrow=37, ncol=3))
colnames(N) = c("Nb", "Nl","Np")
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
  pivot_longer(!year,names_to="species",values_to ="abundance") %>%
  mutate(sqrt_abundance=sqrt(abundance))
Nblp2$species[Nblp2$species == "Nb"] <- "Bromus"
Nblp2$species[Nblp2$species == "Np"] <- "Plantago"
Nblp2$species[Nblp2$species == "Nl"] <- "Layia"

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
  pivot_longer(!year,names_to="species",values_to ="abundance")%>%
  mutate(sqrt_abundance = sqrt(abundance))
Nlp$species[Nlp$species == "Np"] <- "Plantago"
Nlp$species[Nlp$species == "Nl"] <- "Layia"

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

lp <- ggplot(Nlp,aes(year,sqrt_abundance,color=species)) +
  geom_line(size = .8) + xlab("Year") +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),units = "pt"),legend.text = element_markdown(),strip.text.x = element_blank(),panel.spacing = unit(1.5, "lines"),legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
  ylab(expression(Abundance~(m^{"2"})))+
  scale_color_manual(values=c("#0072B2","#009E73"),guide=FALSE)+
  scale_x_continuous(expand = c(0.04,0.04))

sim1 <- ggplot(Nblp,aes(year,sqrt_abundance,color=species)) +
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


sim2 <-  ggplot(Nblp2,aes(year,sqrt_abundance,color=species)) +
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


