library(tidyverse)
library(ggtext)

calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

#all species
dat <- read.csv(paste(datpath, "JR_cover_1mplot.csv", sep = "")) %>%
  filter(treatment == "c") %>%
  filter(species == "PLER" | species == "VUMI" | species == "BRMO" | species == "LAPL") %>%
  group_by(year,species) %>%
  summarize(mean_cov = mean(cover), se_cov = calcSE(cover))

#BRHO, PLER, and LAPL
dat <- read.csv(paste(datpath, "JR_cover_1mplot.csv", sep = "")) %>%
  filter(treatment == "c") %>%
  filter(species == "PLER" | species == "BRMO" | species == "LAPL") %>%
  group_by(year,species) %>%
  summarize(mean_cov = mean(cover), se_cov = calcSE(cover))

theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 13),
              strip.text= element_text(size = 13),
              axis.text = element_text(size = 13))

spp.labs <- c("Bromus", "Layia", "Plantago", "Festuca")
names(spp.labs) <- c("BRMO","LAPL","PLER","VUMI")

#graph with species and treatments
ggplot(dat, aes(year,mean_cov,group=interaction(species,treatment))) + 
  geom_point(aes(color=treatment,shape=treatment),size=2) + 
  geom_line(aes(color=treatment,linetype=treatment),size=1) +
  #geom_errorbar(aes(x=year,y=mean_cov,ymin=mean_cov-se_cov,ymax=mean_cov+se_cov,
  #color = treatment)) +
  scale_linetype_manual(values = c("solid","dashed","longdash"),name="Exclosure treatments", 
                        labels = c("Control","Gopher","Rabbit"))+
  scale_color_manual(values = c("grey40","grey75","black"),name="Exclosure treatments", 
                     labels = c("Control","Gopher","Rabbit")) +
  scale_shape_manual(values = c(16,17,15),name="Exclosure treatments", 
                     labels = c("Control","Gopher","Rabbit")) +
  xlab("Year") + facet_grid(species~.,scale= "free",
                            labeller = labeller(species = spp.labs)) +  
  ylab(expression(Percent~cover~(m^{"2"}))) +
  theme(strip.text = element_text(face = "italic"))

#graph with only species 
ggplot(dat, aes(year,mean_cov,color=species,linetype=species)) + 
  geom_line(size=1) +
  geom_point(size=1) +
  #geom_errorbar(aes(x=year,y=mean_cov,ymin=mean_cov-se_cov,ymax=mean_cov+se_cov))+
  xlab("Year") +  
  ylab(expression(Percent~cover~(m^{"2"}))) +
  theme(legend.text = element_text(face="italic")) +
  scale_linetype_manual(values=c("solid","solid","dashed","solid"),name = "Species", 
                     labels = c("Bromus", "Layia", "Plantago", "Festuca")) +
  scale_color_manual(values=c("gray8","gray47","gray63","gray80"),name = "Species", 
                     labels = c("Bromus", "Layia", "Plantago", "Festuca"))
#without femi
jr <- ggplot(dat, aes(year,mean_cov,color=species)) + 
  geom_line(size=.8) + geom_point(size=1.3) +
  #geom_errorbar(aes(x=year,y=mean_cov,ymin=mean_cov-se_cov,ymax=mean_cov+se_cov))+
  xlab("Year") +  
  ylab(expression(Percent~cover~(m^{"2"}))) +
  theme(legend.text = element_text(face="italic")) +
  scale_color_manual(values=c("grey80","grey50","grey30"),name = "Species", 
                     labels = c("Bromus", "Layia", "Plantago"))
                      
                                                                            