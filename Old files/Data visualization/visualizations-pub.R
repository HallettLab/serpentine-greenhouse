####################################
##### Publication ready graphs #####
####################################

library(tidyverse)
library(ggtext)
library(ggplotify)
library(gridExtra)
library(grid)
library(lattice)
library(gtable)
library(cowplot)

## Read in data
seed_biomass_dat <- read.csv(paste(datpath, "/clean_dat2.csv", sep = ""))
stems_dat <-  read.csv(paste(datpath, "/stems_background.csv", sep = ""))

## SE function
calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

## Set theme
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 12),
              strip.text= element_text(size = 12),
              axis.text = element_text(size = 12))

########################
## Background biomass ##
########################

## Data manipulation 
pler_back <- seed_biomass_dat %>% 
  filter(species == "PLER", background == "PLER")

femi_back <- seed_biomass_dat %>%
  filter(species == "FEMI", background == "FEMI")

lapl_back <- seed_biomass_dat %>%
  filter(species == "LAPL", background == "LAPL")

brho_back <- seed_biomass_dat %>%
  filter(species == "BRHO", background == "BRHO")

femi_pler_back <- full_join(pler_back,femi_back)
brho_lapl_back <- full_join(lapl_back,brho_back)

back_biomass <- full_join(femi_pler_back,brho_lapl_back) %>%
  select(-X,-seeds_in,-seeds_out,-out_in)
  
## Reorder and rename factors
back_biomass$trt_water <- factor(back_biomass$trt_water , levels = c("lo","hi"))
back_biomass$trt_N <- factor(back_biomass$trt_N , levels = c("lo","int","hi"))
back_biomass$background <- factor(back_biomass$background , levels = c("PLER","LAPL","FEMI","BRHO"))

water.labs <- c("Dry", "Wet")
names(water.labs) <- c("lo", "hi")
n.labs <- c("Low","Intermediate","High")
names(n.labs) <- c("lo", "hi")
sp.labs <- c("Bromus hordeaceus", "Festuca microstachys","Layia platyglossa","Plantago erecta")
names(sp.labs) <- c("BRHO", "FEMI","LAPL","PLER")
dens.labs <- c("Low seed addition", "High seed addition")
names(dens.labs) <- c("lo", "hi")

## Data visualization
# V1
ggplot(back_biomass) + geom_boxplot(aes(trt_N,biomass_g,fill=trt_water)) + 
  facet_grid(seed_density~background,
             labeller = labeller(background=sp.labs,seed_density=dens.labs)) +
  ylab("Biomass (g)") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Water treatments", labels = c("Dry","Wet"), 
                    values=c("grey80","grey40")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(strip.text.x = element_text(face = "italic"))

# V2
ggplot(back_biomass) + geom_boxplot(aes(trt_N,biomass_g,fill=background)) + 
  facet_grid(seed_density~trt_water,
             labeller = labeller(trt_water=water.labs,seed_density=dens.labs)) +
  ylab("Biomass (g)") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Background species", labels = c("*Plantago erecta*",
                                                          "*Layia platyglossa*",
                                                          "*Festuca microstachys*",
                                                          "*Bromus hordeaceus*"), 
                    values=c("white","grey85","grey42","grey3")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(legend.text = element_markdown())

#################################
## Background stem/recruitment ##
#################################

## Data manipulation
rec_dat <- stems_dat %>%
  filter(seed_sp != "none") %>%
  mutate(recruitment = (stem_density/seed_added)*100) %>%
  rename(background = seed_sp)

## Reorder and rename factors
rec_dat$trt_water <- factor(rec_dat$trt_water , levels = c("lo","hi"))
rec_dat$trt_N <- factor(rec_dat$trt_N , levels = c("lo","int","hi"))
rec_dat$background <- factor(rec_dat$background , levels = c("PLER","LAPL","FEMI","BRHO"))

## Data visualization
# Recruitment
# V1
ggplot(rec_dat) + geom_boxplot(aes(trt_N,recruitment,fill=trt_water)) + 
  facet_grid(seed_density~background,
             labeller = labeller(background=sp.labs,seed_density=dens.labs)) +
  ylab("Recruitment (%)") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Water treatments", labels = c("Dry","Wet"), 
                    values=c("grey80","grey40")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(strip.text.x = element_text(face = "italic"))

# V2
ggplot(rec_dat) + geom_boxplot(aes(trt_N,recruitment,fill=background)) + 
  facet_grid(seed_density~trt_water,
             labeller = labeller(trt_water=water.labs,seed_density=dens.labs)) +
  ylab("Recruitment(%)") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Background species", labels = c("*Plantago erecta*",
                                                          "*Layia platyglossa*",
                                                          "*Festuca microstachys*",
                                                          "*Bromus hordeaceus*"), 
                    values=c("white","grey85","grey42","grey3")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(legend.text = element_markdown())

# Stem density
# V1
ggplot(rec_dat) + geom_boxplot(aes(trt_N,stem_density,fill=trt_water)) + 
  facet_grid(seed_density~background,
             labeller = labeller(background=sp.labs,seed_density=dens.labs), 
             scale = "free") +
  ylab("Stem density") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Water treatments", labels = c("Dry","Wet"), 
                    values=c("grey80","grey40")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(strip.text.x = element_text(face = "italic"))

# V2
ggplot(rec_dat) + geom_boxplot(aes(trt_N,stem_density,fill=background)) + 
  facet_grid(seed_density~trt_water,
             labeller = labeller(trt_water=water.labs,seed_density=dens.labs),
             scale = "free") +
  ylab("Stem density") + xlab("Nitrogen treatments") + 
  scale_fill_manual(name="Background species", labels = c("*Plantago erecta*",
                                                          "*Layia platyglossa*",
                                                          "*Festuca microstachys*",
                                                          "*Bromus hordeaceus*"), 
                    values=c("white","grey85","grey42","grey3")) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  theme(legend.text = element_markdown())

################################
## Per capita seed production ##
################################

## Data manipulation
seed_dat <- seed_biomass_dat %>%
  filter(seed_density != "none") %>%
  select(-X,-biomass_g) %>%
  group_by(species,background,seed_density,trt_N,trt_water) %>%
  summarize(mean_out_in = mean(out_in,na.rm = T),se_out_in=calcSE(out_in))

seed_dat_none <- seed_biomass_dat %>%
  filter(seed_density == "none") %>%
  select(-X,-biomass_g) %>%
  group_by(species,background,seed_density,trt_N,trt_water) %>%
  summarize(mean_out_in = mean(out_in,na.rm=T),se_out_in=calcSE(out_in))

## Reorder and rename factors
seed_dat$trt_water <- factor(seed_dat$trt_water , levels = c("lo","hi"))
seed_dat$trt_N <- factor(seed_dat$trt_N , levels = c("lo","int","hi"))
seed_dat$background <- factor(seed_dat$background , levels = c("PLER","LAPL","FEMI","BRHO"))
seed_dat$species <- factor(seed_dat$species , levels = c("PLER","LAPL","FEMI","BRHO"))

seed_dat_none$trt_water <- factor(seed_dat_none$trt_water , levels = c("lo","hi"))
seed_dat_none$trt_N <- factor(seed_dat_none$trt_N , levels = c("lo","int","hi"))
seed_dat_none$background <- factor(seed_dat_none$background , levels = c("PLER","LAPL","FEMI","BRHO"))
seed_dat_none$species <- factor(seed_dat_none$species , levels = c("PLER","LAPL","FEMI","BRHO"))

## Data visualization
#For jittering
pd <- position_dodge(0.3)

# Graph of background competition
p <- ggplot(seed_dat, aes(trt_N,mean_out_in,group = interaction(species, trt_water))) +
  geom_point(aes(color = species, shape = trt_water),size=2.5,position=pd) + 
  geom_line(aes(color = species),position=pd) +
  facet_grid(seed_density~background,
             labeller = labeller(background = sp.labs,
                                 seed_density = dens.labs),scale = "free") + 
  ylab("Per capita seed production") +
  scale_shape_manual(values = c(1,16), name="Water treatments", 
                     labels = c("Dry","Wet")) + xlab("Nitrogen treatments") + 
  theme(strip.text.x = element_text(face = "italic"),legend.text = element_markdown()) +
  scale_x_discrete(labels = c("Low","Intermediate","High")) + 
  scale_color_manual(name = "Phytometer species",labels = c("*Bromus*",
                                                             "*Layia*",
                                                             "*Plantago*"),
                       values=c("#D55E00","#0072B2","#009E73")) +
  ggtitle("Background competitor")+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))#+
  #geom_errorbar(aes(x = trt_N, y = mean_out_in, 
                    #ymin = mean_out_in - se_out_in, ymax = mean_out_in + se_out_in, 
                    #color = species),position=pd)

# Graph of no background competition
p2 <- ggplot(seed_dat_none, aes(trt_N,mean_out_in,group = interaction(species, trt_water))) +
  geom_point(aes(color = species, shape = trt_water),size=2.5,position=pd) + 
  geom_line(aes(color = species),position=pd) +
  facet_grid(seed_density~background,scale = "free") + 
  ylab("Per capita seed production") +
  scale_shape_manual(values = c(1,16), name="Water treatments", 
                     labels = c("Dry","Wet")) + xlab("Nitrogen treatments") + 
  theme(strip.text.x = element_blank(),strip.text.y = element_blank(),
        legend.text = element_blank(),legend.position = "none")+#,axis.title.y = element_blank(),axis.title.x = element_blank())+
  scale_x_discrete(labels = c("Low","Intermediate","High")) + 
  ggtitle("No background competitor") +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  scale_color_manual(values=c("#D55E00","#0072B2","#009E73"),name = "Phytometer species",
                     labels = c("*Bromus*","*Layia*","*Plantago*"))

  #geom_errorbar(aes(x = trt_N, y = mean_out_in, 
                    #ymin = mean_out_in - se_out_in, ymax = mean_out_in + se_out_in, 
                    #color = species),position=pd)

## New edition
labelT = "Background species competition"

z <- ggplotGrob(p)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posT <- subset(z$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
height <- z$heights[min(posT$t)]  # height of current top strips

 
z <- gtable_add_rows(z, height, min(posT$t)-1)

# Construct the new strip grobs
stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA)),
  textGrob(labelT, gp = gpar(fontsize = 14, col = "black",fontface="bold"))))

# Position the grobs in the gtable
z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z)

p1 <- as.ggplot(z)

plot_grid(p1,p2,labels=c("A","B"),label_size = 13,
          rel_widths = c(3,1))

ggsave("per_capita_seed.png")

#####################################################
## Per capita seed production other visualizations ##
#####################################################

## Data manipulation
seed_dat <- seed_biomass_dat %>%
  filter(seed_density != "none") %>%
  select(-X,-biomass_g) %>%
  group_by(seed_density,species,background,trt_N,trt_water) %>%
  summarize(mean_out_in = mean(out_in,na.rm=T), se_out_in = calcSE(out_in))

seed_dat_none <- seed_biomass_dat %>%
  filter(seed_density == "none") %>%
  select(-X,-biomass_g) %>% 
  mutate(background = species) %>%
  group_by(species,background,seed_density,trt_N,trt_water) %>%
  summarize(mean_out_in = mean(out_in,na.rm=T),se_out_in=calcSE(out_in))

seed_dat <- full_join(seed_dat,seed_dat_none)

## Reorder and rename factors
seed_dat$trt_water <- factor(seed_dat$trt_water , levels = c("lo","hi"))
seed_dat$trt_N <- factor(seed_dat$trt_N , levels = c("lo","int","hi"))
seed_dat$background <- factor(seed_dat$background , levels = c("PLER","LAPL","FEMI","BRHO"))
seed_dat$species <- factor(seed_dat$species,levels = c("PLER","LAPL","FEMI","BRHO"))
seed_dat$seed_density <- factor(seed_dat$seed_density,levels = c("none","lo","hi"))

water.labs <- c("Dry", "Wet")
names(water.labs) <- c("lo", "hi")
n.labs <- c("Low","Intermediate","High")
names(n.labs) <- c("lo","int","hi")
sp.labs <- c("Plantago erecta","Layia platyglossa","Festuca microstachys","Bromus hordeaceus")
names(sp.labs) <- c("PLER","LAPL","FEMI","BRHO")
dens.labs <- c("None", "Low","High")
names(dens.labs) <- c("none","lo", "hi")

## Data visualization
ggplot(seed_dat, aes(trt_N,mean_out_in,group = interaction(species, trt_water))) +
  geom_point(aes(color = species, shape = trt_water),size=2.5,position=pd) + 
  geom_line(aes(color = species),position=pd) +
  facet_grid(seed_density~background,
             labeller = labeller(background = sp.labs,
                                 seed_density = dens.labs),scale = "free") + 
  ylab("Per capita seed production") +
  scale_shape_manual(values = c(1,16), name="Water treatments", 
                     labels = c("Dry","Wet")) + xlab("N treatments") + 
  theme(strip.text.x = element_text(face = "italic"),legend.text = element_markdown()) +
  #scale_x_discrete(labels = c("Low","Interm","High")) + 
  #scale_color_manual(name = "Background species",labels = c("*Plantago erecta*",
                                                            "*Layia platyglossa*",
                                                            "*Festuca microstachys*",
                                                            "*Bromus hordeaceus*"),
                     values=c("grey80","grey65","grey50","black")) #+
#geom_errorbar(aes(x = trt_N, y = mean_out_in, 
#ymin = mean_out_in - se_out_in, ymax = mean_out_in + se_out_in, 
#color = species),position=pd)

