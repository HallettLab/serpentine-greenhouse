library(tidyverse)
library(ggtext)

calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

#rainfall
ppt <- read.csv(paste(datpath, "JR_rain.csv", sep = "")) %>%
  filter(year != 1982)

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

#graph with species 
jr <- ggplot(dat, aes(year,mean_cov,color=species)) + 
  geom_line(size=.8) + 
  #geom_errorbar(aes(x=year,y=mean_cov,ymin=mean_cov-se_cov,ymax=mean_cov+se_cov))+
  ylab(expression(Percent~cover~(m^{"2"}))) +
  theme(legend.text = element_text(face="italic"),legend.position = "top",axis.title.x = element_blank(),axis.text.x = element_blank())+
  scale_color_manual(values=c("#D55E00","#0072B2","#009E73"),name = "Species", 
                     labels = c("Bromus", "Layia", "Plantago"))+
    scale_x_continuous(expand = c(0.04, 0.04)) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm"))


p <- ggplot(ppt, aes(year,growing_season_ppt)) + 
  geom_line(size=.8) +
  xlab("Year") +  
  ylab("Precipitation (mm)") +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm"))


gjr <- ggplotGrob(jr)
gp <- ggplotGrob(p)
maxWidth = grid::unit.pmax(gjr$widths[2:5], gp$widths[2:5])
gjr$widths[2:5] <- as.list(maxWidth)
gp$widths[2:5] <- as.list(maxWidth)
grid.arrange(gjr,gp,ncol=1,left = textGrob(c("a)","b)"), x =c(2.7,2.7), 
                                    y = c(.92,.51), gp = gpar(fontface = "bold", fontsize = 13)))
