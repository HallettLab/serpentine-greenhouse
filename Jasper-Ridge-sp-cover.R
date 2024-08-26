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
dat <- read.csv(paste(datpath, "JR_cover_1mplot1983-2019.csv",sep="")) %>%
  mutate_at(4,as.character) %>%
  filter(treatment=="c") %>%
  filter(species=="BRMO"|species=="LAPL"|species =="LOMU"|species=="PLER") %>%
  mutate(species = case_when(species %in% "BRMO" ~ "Bromus",
                             species %in% "LAPL" ~ "Layia",
                             species %in% "PLER" ~ "Plantago"))

jrdat <-  dat %>%
  group_by(year,species) %>%
  summarize(mean_cov = mean(cover),med_cov = median(cover))
CI <- dat %>%
  group_by(year,species) %>%
  summarize(lower=quantile(cover,probs = 0.25),upper=quantile(cover,probs=0.75))
jrdatCI <- left_join(jrdat,CI)

#graph with species 
#without lomu for ms
jrms <- jrdatCI %>%
  filter(year %in% 1983:2019)

jrmsfig <- ggplot(jrms,aes(year,med_cov,color=species)) +
  geom_ribbon(aes(x= year, ymin=lower,ymax=upper, fill=species),colour=NA, alpha = .2) +
  geom_line(linewidth=.8) +
  ylab("Percent cover") +
  scale_color_manual(values=c("#D55E00","#0072B2","#48B99B"),name = "Species", 
                     labels = c("*Bromus*", "*Layia*","*Plantago*"))+
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                    values=c("#D55E00","#0072B2","#48B99B")) +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  xlab("Year") +
  theme(legend.position = "none",legend.text = ggtext::element_markdown()) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm")) 

#rainfall
ppt <- read.csv(paste(datpath, "JR_rain.csv", sep = "")) %>%
  filter(year %in% c(1983:2019)) 

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

#2019
gjr <- ggplotGrob(jrmsfig)
gp <- ggplotGrob(p2019)
maxWidth = grid::unit.pmax(gjr$widths[2:5], gp$widths[2:5])
gjr$widths[2:5] <- as.list(maxWidth)
gp$widths[2:5] <- as.list(maxWidth)
grid.arrange(gjr,gp,ncol=1,left = textGrob(c("a)","b)"), x =c(2.7,2.7), 
                                           y = c(.92,.51), gp = gpar(fontface = "bold", fontsize = 16)))

pdf("jr-timeseries-msfig.pdf")
plot_grid(jrmsfig,p2019,ncol=1,align="v",rel_heights = c(.2,.09),hjust = -3,vjust=0)
dev.off()

##legend if needed
legmsfig <- ggplot(jrms,aes(year,med_cov,color=species)) + 
  geom_line(size=.8) + 
  geom_ribbon(aes(x= year, ymin=lower,ymax=upper, fill=species), alpha = .2) +
  ylab("Percent cover") +
  scale_color_manual(values=c("#D55E00","#0072B2","#48B99B"),name = "Species", 
                     labels = c("*Bromus*", "*Layia*","*Plantago*"))+
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Plantago*"),
                    values=c("#D55E00","#0072B2","#48B99B")) +
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

legendms <- as.ggplot(g_legend(legmsfig))




#####################################################
#with lomu

#skip to line 78 to read in data by code produced below
#1983-2019 1m^s plot cover
#convert 2021-23 cover to 1m^2 plot 
dat_upd <- read.csv(paste(datpath, "JR_cover2021-23.csv", sep = ""))

cover_upd <- dat_upd %>%
  separate(1,into=c("treatment","trtrep"),sep=1) %>%
  mutate(treatment = "c") %>%
  filter(treatment == "c") %>%
  rename(year=Year,species=Species,cover=Data,subplot=Quad) %>%
  select(-8) %>%
  filter(species == "Plantago" | species == "Bromus" | species == "Layia" | species == "Lolium") %>%
  mutate_at(6,as.numeric) %>%
  mutate_at(4,as.character) 

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
JRm2cover <- merge(cover_upd, key) %>%
  tbl_df() %>%
  group_by(year, species, treatment, trtrep, plot) %>%
  summarize(cover=mean(cover)) %>%
  mutate(uniqueID=paste(treatment, trtrep, sep="_")) %>%
  #removing the middle block for independence
  filter(plot!="block1_middle", plot!="block2_middle") %>%
  mutate(treatment="c") 

cover_complete <- rbind(dat,JRm2cover)
jrlomu <- ggplot(jrdatCI,aes(year,med_cov,color=species)) + 
  geom_line(size=.8) + 
  geom_ribbon(aes(x= year, ymin=lower,ymax=upper, fill=species), alpha = .2) +
  ylab("Percent cover") +
  xlab("Year") +
  scale_color_manual(values=c("#D55E00","#0072B2","black","#48B99B"),name = "Species", 
                     labels = c("*Bromus*", "*Layia*","*Lolium*","*Plantago*"))+
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Lolium*","*Plantago*"),
                    values=c("#D55E00","#0072B2","black","#48B99B")) +
  theme(legend.position = "none",legend.text = ggtext::element_markdown()) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm"))



leglomu <- ggplot(jrdatCI,aes(year,med_cov,color=species)) + 
  geom_line(size=.8) + 
  geom_ribbon(aes(x= year, ymin=lower,ymax=upper, fill=species), alpha = .2) +
  ylab("Percent cover") +
  xlab("Year") +
  scale_color_manual(values=c("#D55E00","#0072B2","black","#48B99B"),name = "Species", 
                     labels = c("*Bromus*", "*Layia*","*Lolium*","*Plantago*"))+
  scale_fill_manual(name = "Species",labels = c("*Bromus*","*Layia*","*Lolium*","*Plantago*"),
                    values=c("#D55E00","#0072B2","black","#48B99B")) +
  theme(legend.position = "top",legend.text = ggtext::element_markdown()) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm"))

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

legendlomu <- as.ggplot(g_legend(leglomu))


p2023 <- ggplot(ppt2023, aes(x=year)) + 
  geom_line(aes(y=growing_season_ppt),data=ppt[1:5,],colour="black",size=.8) +
  geom_line(aes(y=growing_season_ppt),data=ppt[5:9,],colour="#c54e49",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[9:25,],colour="black",size=.8) +
  geom_line(aes(y=growing_season_ppt),data=ppt[25:27,],colour="#c54e49",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[27:30,],colour="black",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[30:33,],colour="#c54e49",size=.8)+
  geom_line(aes(y=growing_season_ppt),data=ppt[33:41,],colour="black",size=.8)+
  xlab("Year") + 
  #ylab(expression(atop("Rainfall",paste("(mm)"))))+
  ylab("Rainfall (mm)") +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  theme(plot.margin = unit(c(.01,.8,.5,.5), "cm")) +
  geom_hline(yintercept=565, linetype="dashed")+
  annotate("pointrange", x =1983, y =1250.442, 
           ymin = 1250.442, ymax = 1250.442,
           colour = "#4b7ea3")+
  annotate("pointrange", x =1998, y =1028.446, 
           ymin = 1028.446, ymax = 1028.446,
           colour = "#4b7ea3")+
  annotate("pointrange", x=2017,y=859.028,ymin = 859.028, ymax = 859.028,
           colour = "#4b7ea3")+
  annotate("pointrange", x=2023,y=1537.54,ymin = 1537.54, ymax = 1537.54,
           colour = "#4b7ea3") 

#2023
gjrlomu <- ggplotGrob(jrlomu)
gp <- ggplotGrob(p2023)
maxWidth = grid::unit.pmax(gjr$widths[2:5], gp$widths[2:5])
gjr$widths[2:5] <- as.list(maxWidth)
gp$widths[2:5] <- as.list(maxWidth)
grid.arrange(gjrlomu,gp,ncol=1,left = textGrob(c("a)","b)"), x =c(2.7,2.7), 
                                               y = c(.92,.51), gp = gpar(fontface = "bold", fontsize = 16)))

pdf("jr-timeseries-lomu.pdf", width = , height = )
plot_grid(jrmsfig,p2023,ncol=1,align="v",rel_heights = c(.03,.2,.05),hjust = -3,vjust=0)
dev.off()

