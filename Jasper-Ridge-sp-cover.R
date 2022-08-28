library(tidyverse)
library(ggtext)
library(gridExtra)
library(grid)

calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 16),
              axis.text = element_text(size = 16))

#BRHO, PLER, and LAPL
dat <- read.csv(paste(datpath, "JR_cover_1mplot.csv", sep = "")) %>%
  filter(treatment == "c") %>%
  filter(species == "PLER" | species == "BRMO" | species == "LAPL") %>%
group_by(year,species) %>%
  summarize(mean_cov = mean(cover), se_cov = calcSE(cover))

spp.labs <- c("Bromus", "Layia", "Plantago")
names(spp.labs) <- c("BRMO","LAPL","PLER")

#graph with species 
jr <- ggplot(dat,aes(year,mean_cov,color=species)) + 
  geom_line(size=.8) + 
  geom_errorbar(aes(x=year,y=mean_cov,ymin=mean_cov-se_cov,ymax=mean_cov+se_cov))+
  ylab(expression(Percent~cover~(m^{"2"}))) +
  scale_color_manual(values=c("#D55E00","#0072B2","#009E73"),name = "Species", 
                     labels = c("*Bromus*", "*Layia*", "*Plantago*"))+
    scale_x_continuous(expand = c(0.04, 0.04)) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm"), legend.position = "top")+
  xlab("Year") +
  theme(legend.text = ggtext::element_markdown()) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())


#rainfall
ppt <- read.csv(paste(datpath, "JR_rain.csv", sep = "")) %>%
  filter(year != 1982)

x <- expression(paste("El Ni", tilde(n), "o"))

p <- ggplot(ppt, aes(x=year)) + 
  geom_line(aes(y=growing_season_ppt),data=ppt[1:5,],colour="black",size=.8) +
  geom_line(aes(y=growing_season_ppt),data=ppt[5:9,],colour="#de2d26",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[9:25,],colour="black",size=.8) +
  geom_line(aes(y=growing_season_ppt),data=ppt[25:27,],colour="#de2d26",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[27:30,],colour="black",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[30:33,],colour="#de2d26",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[33:37,],colour="black",size=.8)+
  xlab("Year") +  
  ylab("Precipitation (mm)") +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm")) +
  geom_hline(yintercept=565, linetype="dashed")+
  annotate("pointrange", x =1983, y =1250.442, 
           ymin = 1250.442, ymax = 1250.442,
           colour = "#3182bd")+
  annotate("text",x=1983.5,y=1292,label=x,colour="#3182bd")+
  annotate("pointrange", x =1998, y =1028.446, 
           ymin = 1028.446, ymax = 1028.446,
           colour = "#3182bd")+
  annotate("text", x=1998,y=1072,label=x,colour="#3182bd") +
  annotate("pointrange", x=2017,y=859.028,ymin = 859.028, ymax = 859.028,
           colour = "#3182bd")+
  annotate("text", x=2017,y=902,label=x,colour="#3182bd") +
  annotate("text",x=2014,y=205,label="Drought",colour="#de2d26")+
  annotate("text",x=2008.5,y=268,label="Drought",colour="#de2d26")+
  annotate("text",x=1989,y=275,label="Drought",colour="#de2d26")

gjr <- ggplotGrob(jr)
gp <- ggplotGrob(p)
maxWidth = grid::unit.pmax(gjr$widths[2:5], gp$widths[2:5])
gjr$widths[2:5] <- as.list(maxWidth)
gp$widths[2:5] <- as.list(maxWidth)
grid.arrange(gjr,gp,ncol=1,left = textGrob(c("a)","b)"), x =c(2.7,2.7), 
                                    y = c(.92,.51), gp = gpar(fontface = "bold", fontsize = 16)))
