x <- expression(paste("El Ni", tilde(n), "o"))

p <- ggplot(ppt, aes(year,growing_season_ppt)) + 
  geom_line(size=.8) +
  xlab("Year") +  
  ylab("Precipitation (mm)") +
  scale_x_continuous(expand = c(0.04, 0.04)) +
  theme(plot.margin = unit(c(.5,.8,.5,.5), "cm")) +
  geom_hline(yintercept=565, linetype="dashed")+
  annotate("text",x=1983,y=1275,label=x)+
  annotate("text", x=1998,y=1055,label=x) +
  annotate("text", x=2017,y=885,label=x) +
  annotate("label",x=2014,y=200,label="Prolonged drought")+
  annotate("label",x=2008,y=250,label="Prolonged\ndrought")+
  annotate("label",x=1989,y=275,label="Prolonged drought")

