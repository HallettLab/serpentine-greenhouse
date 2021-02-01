seed_dat_nodens <- seed_dat %>%
  filter(individual == "LAPL") %>%
  group_by(individual,background_comp,trt_water,trt_N) %>%
  summarize(mean_out_in = mean(out_in, na.rm = T),se_out_in = calcSE(out_in))

seed_dat_nodens$trt_water <- factor(seed_dat_nodens$trt_water , levels = c("lo","hi"))
seed_dat_nodens$trt_N <- factor(seed_dat_nodens$trt_N , levels = c("lo","int","hi"))

water.labs <- c("Low water", "High water")
names(water.labs) <- c("lo", "hi")

ggplot(seed_dat_nodens, aes(trt_N,mean_out_in, group = interaction(individual, trt_water))) +
  theme_bw() + geom_point(aes(color = background_comp)) + geom_line(aes(color = background_comp)) + ylab("Per capita seed production") +
  facet_grid(individual~trt_water, labeller = labeller(trt_water = water.labs), scale = "free")
  
ggplot(seed_dat_nodens, aes(trt_N,mean_out_in)) + geom_point(aes(color = background_comp)) +
  facet_wrap(~trt_water)

seed_dat <- read.csv(paste(datpath, "/clean_seed_dat.csv", sep = ""))

seed_dat <- seed_dat %>%
  group_by(seed_density,individual,background_comp,trt_water,trt_N) %>%
  summarize(mean_out_in = mean(out_in, na.rm = T),se_out_in = calcSE(out_in))

seed_dat$trt_water <- factor(seed_dat$trt_water , levels = c("lo","hi"))
seed_dat$trt_N <- factor(seed_dat$trt_N , levels = c("lo","int","hi"))

ggplot(subset(seed_dat, !seed_density %in% c("none","lo") & individual %in% "LAPL"), aes(trt_N,mean_out_in,group = interaction(background_comp, seed_density))) + 
  geom_point(aes(color = background_comp, shape = seed_density),size=2.5) + geom_line(aes(color = background_comp)) +
  facet_grid(individual~trt_water,scale = "free") + ylab("Per capita seed production") +
  theme_bw() + ylim(0,20)


ggplot(subset(seed_dat, !seed_density %in% c("none","hi") & individual %in% "PLER" & !background_comp %in% "PLER"), aes(trt_N,mean_out_in,group = interaction(background_comp, seed_density))) + 
  geom_point(aes(color = background_comp, shape = seed_density),size=2.5) + geom_line(aes(color = background_comp)) +
  facet_grid(individual~trt_water,scale = "free") + ylab("Per capita seed production") + 
  theme_bw() + ylim(0,20)


ggplot(subset(seed_dat, seed_density %in% "none" & individual %in% "LAPL"), aes(trt_N,mean_out_in,group = interaction(background_comp, seed_density))) + 
  geom_point(aes(color = background_comp, shape = seed_density),size=2.5) + geom_line(aes(color = background_comp)) +
  facet_grid(individual~trt_water,scale = "free") + ylab("Per capita seed production") + theme_bw()

ggsave("none_dens.png")

seed_pler <- seed_dat %>%
  filter(individual == "PLER", seed_density == "hi")  %>%
  filter(background_comp == "BRHO" | background_comp == "LAPL")
  
ggplot(seed_pler, aes(trt_N,mean_out_in,group = interaction(background_comp, seed_density))) + 
  geom_point(aes(color = background_comp, shape = seed_density),size=2.5) + geom_line(aes(color = background_comp)) +
  facet_grid(individual~trt_water,scale = "free") + ylab("Per capita seed production") 
