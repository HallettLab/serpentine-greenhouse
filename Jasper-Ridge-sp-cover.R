library(tidyverse)
library(ggtext)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplotify)

theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 16),
              axis.text = element_text(size = 16))

#####################################
######Code for updated JR data#######
#####################################
#skip to line 52 to read in data produced by code below
#BRHO, PLER, and LAPL 1983-2019
dat <- read.csv(paste(datpath, "JR_cover_1mplot1983-2019.csv", sep = "")) %>%
  filter(treatment == "c") %>%
  filter(species == "PLER" | species == "BRMO" | species == "LAPL") %>%
  mutate(species=case_when(species %in% "BRMO"~"Bromus",
                           species %in%"LAPL" ~ "Layia",
                           species %in% "PLER" ~"Plantago")) %>%
  mutate_at(4,as.character)
#BRHO, PLER, and LAPL 2021-2021
cover <- read.csv(paste(datpath,"JR_cover2021-23.csv",sep=""))

cov <- cover %>%
  separate(1,into=c("treatment","trtrep"),sep=-1) %>%
  rename(year=Year,species=Species,subplot=Quad,cover=Data) %>%
  filter(species %in% c("Bromus","Layia","Plantago")) %>%
  mutate_at(4,as.character) %>%
  mutate_at(6,as.numeric) %>%
  select(-8)

# create a key to aggregate 0.5 m x 0.5 m subplots into 1 m x 1 m plots
top_block1 <- cbind(c(1,2,7,8), rep("block1_top", 4))
top_block2 <- cbind(c(13, 14,19, 20), rep("block2_top", 4))

middle_block1 <- cbind(c(3,4,9,10), rep("block1_middle", 4))
middle_block2 <- cbind(c(15, 16, 21, 22), rep("block2_middle", 4))

bottom_block1 <- cbind(c(5, 6, 11, 12), rep("block1_bottom", 4))
bottom_block2 <- cbind(c(17,18, 23, 24 ), rep("block2_bottom", 4))

key <-as.data.frame(rbind(top_block1, top_block2, 
                          middle_block1, middle_block2, 
                          bottom_block1, bottom_block2))
names(key)=c("subplot", "plot")

key <- key %>%
  tbl_df() %>%
  mutate(subplot=as.numeric(as.character(subplot)), 
         plot=as.character(plot))

rm(top_block1, top_block2, middle_block1, middle_block2, bottom_block1, bottom_block2)


# merge JRcover and key; aggregate at the 1 m x 1 m plot scale
JRm2cover <- merge(cov, key) %>%
  tbl_df() %>%
  group_by(year, species, treatment, trtrep, plot) %>%
  summarize(cover=mean(cover)) %>%
  mutate(uniqueID=paste(treatment, trtrep, plot, sep="_")) %>%
  #removing the middle block for independence
  filter(plot!="block1_middle", plot!="block2_middle") %>%
  mutate(treatment="c") %>%
  rename()

comp_cover <- rbind(dat,JRm2cover)

#write.csv(comp_cover,"JR_cover_1mplot1983-2023.csv")
#################################################################

comp_cover <- read.csv(paste(datpath,"JR_cover_1mplot1983-2023.csv",sep=""))

jrdat <-  comp_cover %>%
  group_by(year,species) %>%
  summarize(mean_cov = mean(cover))
lower <- comp_cover %>%
  group_by(year,species) %>%
  summarize(lower=quantile(cover,probs = 0.25))
upper <- comp_cover %>%
  group_by(year,species) %>%
  summarize(upper=quantile(cover,probs=0.75))
CI <- left_join(lower,upper)
jrdatCI <- left_join(jrdat,CI)


#graph with species 

jr <- ggplot(jrdatCI,aes(year,mean_cov,color=species)) + 
  geom_line(size=.8) + 
  geom_ribbon(aes(x= year, ymin=lower,ymax=upper, fill=species), alpha = .2) +
  ylab("Percent cover") +
  scale_color_manual(values=c("#D55E00","#0072B2","#009E73"),name = "Species", 
                     labels = c("*Bromus*", "*Layia*", "*Plantago*"))+
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                    values=c("#D55E00","#0072B2","#009E73")) +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  xlab("Year") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm"))#, legend.position = "top")


leg <- ggplot(jrdatCI,aes(year,mean_cov,color=species)) + 
  geom_line(size=.8) + 
  geom_ribbon(aes(x= year, ymin=lower,ymax=upper, fill=species), alpha = .2) +
  ylab("Percent cover") +
  scale_color_manual(values=c("#D55E00","#0072B2","#009E73"),name = "Species", 
                     labels = c("*Bromus*", "*Layia*", "*Plantago*"))+
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                    values=c("#D55E00","#0072B2","#009E73")) +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  xlab("Year") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm")) +
  theme(legend.text = ggtext::element_markdown())
  

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

legend <- as.ggplot(g_legend(leg))



#rainfall
ppt <- read.csv(paste(datpath, "JR_rain.csv", sep = "")) %>%
  filter(year != 1982)

p2 <- ggplot(ppt, aes(x=year)) + 
  geom_line(aes(y=growing_season_ppt),data=ppt[1:5,],colour="black",size=.8) +
  geom_line(aes(y=growing_season_ppt),data=ppt[5:9,],colour="#de2d26",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[9:25,],colour="black",size=.8) +
  geom_line(aes(y=growing_season_ppt),data=ppt[25:27,],colour="#de2d26",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[27:30,],colour="black",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[30:33,],colour="#de2d26",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[33:41,],colour="black",size=.8)+
  xlab("Year") + 
  #ylab(expression(atop("Rainfall",paste("(mm)"))))+
  ylab("Rainfall (mm)") +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  theme(plot.margin = unit(c(.01,.8,.5,.5), "cm")) +
  geom_hline(yintercept=565, linetype="dashed")+
  annotate("pointrange", x =1983, y =1250.442, 
           ymin = 1250.442, ymax = 1250.442,
           colour = "#3182bd")+
  annotate("pointrange", x =1998, y =1028.446, 
           ymin = 1028.446, ymax = 1028.446,
           colour = "#3182bd")+
  annotate("pointrange", x=2017,y=859.028,ymin = 859.028, ymax = 859.028,
           colour = "#3182bd")+
  annotate("pointrange", x=2023,y=1537.54,ymin = 1537.54, ymax = 1537.54,
           colour = "#3182bd") +
  scale_y_continuous(breaks = seq(0,1600,by=565))

gjr <- ggplotGrob(jr)
gp <- ggplotGrob(p2)
maxWidth = grid::unit.pmax(gjr$widths[2:5], gp$widths[2:5])
gjr$widths[2:5] <- as.list(maxWidth)
gp$widths[2:5] <- as.list(maxWidth)
grid.arrange(gjr,gp,ncol=1,left = textGrob(c("a)","b)"), x =c(2.7,2.7), 
                                    y = c(.92,.51), gp = gpar(fontface = "bold", fontsize = 16)))

plot_grid(legend,jr,p2,ncol=1,align="v",rel_heights = c(.05,.6,.15),labels = c("","a)","b)"),hjust = -3,vjust=0)


