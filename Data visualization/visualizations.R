##########################
######Visualization#######
##########################

library(tidyverse)
library(ggplotify)
library(gridExtra)

## Read in data
seed_dat <- read.csv(paste(datpath, "/clean_seed_dat.csv", sep = ""))

## SE function
calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

## Data manipulation
seed_dat <- seed_dat %>%
  group_by(block,seed_density,individual,background_comp,trt_water,trt_N) %>%
  summarize(out_in = out_in)

dat_exc_none <- dat %>%
  filter(background_comp != "none")

dat_none <- dat %>%
  filter(background_comp == "none")

## Data visualization

# By species graphs
dat_exc_none$trt_water <- factor(dat_exc_none$trt_water , levels = c("lo","hi"))
dat_exc_none$trt_N <- factor(dat_none$trt_N , levels = c("lo","int","hi"))
dat_none$trt_water <- factor(dat_none$trt_water , levels = c("lo","hi"))
dat_none$trt_N <- factor(dat_none$trt_N , levels = c("lo","int","hi"))

density.labs <- c("Low", "High")
names(density.labs) <- c("lo", "hi")

spp.labs <- c("B. hordeaceus", "F. microstachys", "L. platyglossa", "P. erecta")
names(spp.labs) <- c("BRHO","FEMI","LAPL","PLER")


ggplot(subset(dat_exc_none, individual %in% "BRHO")) + geom_boxplot(aes(trt_N,out_in,fill=trt_water)) +
  theme_bw() + facet_grid(seed_density~background_comp, labeller = labeller(seed_density = density.labs, background_comp = spp.labs)) + ylab("*B. hordeaceus* per capita seed production") +
  xlab(expression(Nitrogen~treatments~(kg~N~ha^{"-1"}~year^{"-1"}))) + theme(strip.text.x = element_text(face = "italic")) +
  theme(axis.title.y = ggtext::element_markdown()) + scale_x_discrete(labels = c("0","11","55")) + scale_fill_manual(name="Water treatments", labels = c("Low","High"), values=c("azure3","azure4")) +
  theme(legend.position="bottom") 

ggplot(subset(dat_exc_none, individual %in% "FEMI")) + geom_boxplot(aes(trt_N,out_in,fill=trt_water)) +
  theme_bw() + facet_grid(seed_density~background_comp, labeller = labeller(seed_density = density.labs, background_comp = spp.labs)) + ylab("*F. microstachys* per capita seed production") +
  xlab(expression(Nitrogen~treatments~(kg~N~ha^{"-1"}~year^{"-1"}))) + theme(strip.text.x = element_text(face = "italic")) +
  theme(axis.title.y = ggtext::element_markdown()) + scale_x_discrete(labels = c("0","11","55")) + scale_fill_manual(name="Water treatments", labels = c("Low","High"), values=c("azure3","azure4")) +
  theme(legend.position="bottom") 

ggplot(subset(dat_exc_none, individual %in% "LAPL")) + geom_boxplot(aes(trt_N,out_in,fill=trt_water)) +
  theme_bw() + facet_grid(seed_density~background_comp, labeller = labeller(seed_density = density.labs, background_comp = spp.labs)) + ylab("*L. platyglossa* per capita seed production") +
  xlab(expression(Nitrogen~treatments~(kg~N~ha^{"-1"}~year^{"-1"}))) + theme(strip.text.x = element_text(face = "italic")) +
  theme(axis.title.y = ggtext::element_markdown()) + scale_x_discrete(labels = c("0","11","55")) + scale_fill_manual(name="Water treatments", labels = c("Low","High"), values=c("azure3","azure4")) +
  theme(legend.position="bottom") 

ggplot(subset(dat_exc_none, individual %in% "PLER")) + geom_boxplot(aes(trt_N,out_in,fill=trt_water)) +
  theme_bw() + facet_grid(seed_density~background_comp, labeller = labeller(seed_density = density.labs, background_comp = spp.labs)) + ylab("*P. erecta* per capita seed production") +
  xlab(expression(Nitrogen~treatments~(kg~N~ha^{"-1"}~year^{"-1"}))) + theme(strip.text.x = element_text(face = "italic")) +
  theme(axis.title.y = ggtext::element_markdown()) + scale_x_discrete(labels = c("0","11","55")) + scale_fill_manual(name="Water treatments", labels = c("Low","High"), values=c("azure3","azure4")) +
  theme(legend.position="bottom") 
  
z <- ggplotGrob(p)
# New strip at the top
z <- ggplotGrob(p)
z <- gtable_add_rows(z, z$height[7], pos = 6)
# New strip spans columns 5 to 9
z <- gtable_add_grob(z, 
                     list(rectGrob(gp = gpar(col = NA, size = 0.5)),
                          textGrob("Background species competition", gp = gpar(cex = 1, fontface = 'bold', col = "black"))), 
                     t=7, l=5, b=7, r=11, name = c("a", "b"))
# Add small gap between strips - below row 6
z <- gtable_add_rows(z, unit(2/10, "line"), 7)


# New strip on the right
text = "Seeding density"
size = 13
face = "bold"

# Get the positions of the right strips in the layout: t = top, l = left, ...
strip <-c(subset(z$layout, grepl("strip-r", z$layout$name), select = t:r))

# Text grob
text.grob = textGrob(text, rot = -90,
                     gp = gpar(fontsize = size, fontface = face))

# New column to the right of current strip
# Adjusts its width to text size
width = unit(2, "grobwidth", text.grob) + unit(1, "lines")
z <- gtable_add_cols(z, width, max(strip$r))  

z <- gtable_add_grob(z, text.grob, 
                     t = 2, l = max(strip$r) + 1, b = 14)

# Draw it
grid.newpage()
grid.draw(z)

p1 <- as.ggplot(z)
p2 <- ggplot(subset(dat_none, individual %in% "PLER")) + geom_boxplot(aes(trt_N,out_in,fill=trt_water)) +
  theme_bw() + scale_x_discrete(labels=c("0","11","55")) + scale_fill_manual(values=c("azure3","azure4")) +
  ggtitle("No background competition") + theme(plot.title = element_text(hjust = 0.5,face="bold"),legend.position = 'none',
                                               axis.title = element_blank())

#xlab(expression(Nitrogen~treatments~(kg~N~ha^{"-1"}~year^{"-1"})))
#ylab("*B. hordeaceus* per capita seed production") + theme(axis.title.y = ggtext::element_markdown())
#name="Water treatments", labels = c("Low","High"),

plot_grid(p1,p2,labels=c("A","B"),label_size = 13,
          rel_widths = c(3,1))

## All species graphs
seed_dat <- seed_dat %>%
  group_by(seed_density,individual,background_comp,trt_water,trt_N) %>%
  summarize(mean_out_in = mean(out_in, na.rm = T),se_out_in = calcSE(out_in))

dat_exc_none <- dat_exc_none %>%
  group_by(seed_density,individual,background_comp,trt_water,trt_N) %>%
  summarize(mean_out_in = mean(out_in, na.rm = T),se_out_in = calcSE(out_in)) %>%
  filter(background_comp != "none") 
  
dat_none <- dat_none %>%
  group_by(seed_density,individual,background_comp,trt_water,trt_N) %>%
  summarize(mean_out_in = mean(out_in, na.rm = T),se_out_in = calcSE(out_in))

ggplot(dat_exc_none, aes(trt_N,mean_out_in,group = interaction(individual, trt_water))) +
  theme_bw() + geom_point(aes(color = individual, shape = trt_water),size=2.5) + geom_line(aes(color = individual)) +
  facet_grid(seed_density~background_comp, labeller = labeller(background_comp = spp.labs,seed_density = density.labs),scale = "free") + ylab("Per capita seed production") +
  scale_shape_manual(values = c(1,16), name="Water treatments", labels = c("Low","High")) + xlab(expression(Nitrogen~treatments~(kg~N~ha^{"-1"}~year^{"-1"}))) + theme(strip.text.x = element_text(face = "italic")) +
  scale_x_discrete(labels = c("0","11","55")) + scale_color_discrete(name = "Species phytometer") +
  geom_errorbar(aes(x = trt_N, y = mean_out_in, ymin = mean_out_in - se_out_in, ymax = mean_out_in + se_out_in, 
                    color = individual), width = 0.2)


ggplot(dat_none, aes(trt_N,mean_out_in,group = interaction(individual, trt_water))) +
  theme_bw() + geom_point(aes(color = individual, shape = trt_water),size=2.5) + geom_line(aes(color = individual)) + ylab("Per capita seed production") +
  scale_shape_manual(values = c(1,16), name="Water treatments", labels = c("Low","High")) + xlab(expression(Nitrogen~treatments~(kg~N~ha^{"-1"}~year^{"-1"}))) + theme(strip.text.x = element_text(face = "italic")) +
  scale_x_discrete(labels = c("0","11","55")) + scale_color_discrete(name = "Species phytometer") + ggtitle("No background competition") +
  theme(plot.title = element_text(hjust = 0.5,face="bold"),legend.position = 'none',axis.title = element_blank()) +
  geom_errorbar(aes(x = trt_N, y = mean_out_in, ymin = mean_out_in - se_out_in, ymax = mean_out_in + se_out_in, 
                    color = individual), width = 0.2)




