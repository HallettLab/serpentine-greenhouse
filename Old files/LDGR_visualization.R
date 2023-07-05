## skip to line 45 to read-in data created from the code below

mean_ldgr <- data.frame(rbind(brho_lapl_hi_hi_mean,brho_lapl_hi_int_mean,brho_lapl_hi_lo_mean,
                   brho_lapl_lo_hi_mean,brho_lapl_lo_int_mean,brho_lapl_lo_lo_mean,
                   brho_pler_hi_hi_mean,brho_pler_hi_int_mean,brho_pler_hi_lo_mean,
                   brho_pler_lo_hi_mean,brho_pler_lo_int_mean,brho_pler_lo_lo_mean,
                   pler_lapl_hi_hi_mean,pler_lapl_hi_int_mean,pler_lapl_hi_lo_mean,
                   pler_lapl_lo_hi_mean,pler_lapl_lo_int_mean,pler_lapl_lo_lo_mean,
                   pler_brho_hi_hi_mean,pler_brho_hi_int_mean,pler_brho_hi_lo_mean,
                   pler_brho_lo_hi_mean,pler_brho_lo_int_mean,pler_brho_lo_lo_mean,
                   lapl_brho_hi_hi_mean,lapl_brho_hi_int_mean,lapl_brho_hi_lo_mean,
                   lapl_brho_lo_hi_mean,lapl_brho_lo_int_mean,lapl_brho_lo_lo_mean,
                   lapl_pler_hi_hi_mean,lapl_pler_hi_int_mean,lapl_pler_hi_lo_mean,
                   lapl_pler_lo_hi_mean,lapl_pler_lo_int_mean,lapl_pler_lo_lo_mean))
colnames(mean_ldgr) <- ("mean_ldgr")


sd_ldgr <- data.frame(rbind(brho_lapl_hi_hi_sd,brho_lapl_hi_int_sd,brho_lapl_hi_lo_sd,
                            brho_lapl_lo_hi_sd,brho_lapl_lo_int_sd,brho_lapl_lo_lo_sd,
                            brho_pler_hi_hi_sd,brho_pler_hi_int_sd,brho_pler_hi_lo_sd,
                            brho_pler_lo_hi_sd,brho_pler_lo_int_sd,brho_pler_lo_lo_sd,
                            pler_lapl_hi_hi_sd,pler_lapl_hi_int_sd,pler_lapl_hi_lo_sd,
                            pler_lapl_lo_hi_sd,pler_lapl_lo_int_sd,pler_lapl_lo_lo_sd,
                            pler_brho_hi_hi_sd,pler_brho_hi_int_sd,pler_brho_hi_lo_sd,
                            pler_brho_lo_hi_sd,pler_brho_lo_int_sd,pler_brho_lo_lo_sd,
                            lapl_brho_hi_hi_sd,lapl_brho_hi_int_sd,lapl_brho_hi_lo_sd,
                            lapl_brho_lo_hi_sd,lapl_brho_lo_int_sd,lapl_brho_lo_lo_sd,
                            lapl_pler_hi_hi_sd,lapl_pler_hi_int_sd,lapl_pler_hi_lo_sd,
                            lapl_pler_lo_hi_sd,lapl_pler_lo_int_sd,lapl_pler_lo_lo_sd))
colnames(sd_ldgr) <- ("sd_ldgr")

mean_sd_ldgr <- (cbind(mean_ldgr,sd_ldgr))

write.csv(mean_sd_ldgr,"mean_sd_ldgr.csv")

mean_sd_ldgr <- read.csv(paste(datpath, "mean_sd_ldgr.csv", sep = ""))
colnames(mean_sd_ldgr) <- c("treatment","mean","stdv")

library(tidyverse)

mean_sd_ldgr <- mean_sd_ldgr %>%
  separate(treatment,into=c("invader","resident","water","N"),sep="_")%>%
  unite("invader_resident",1:2)

write.csv(mean_sd_ldgr,"mean_sd_ldgr.csv")

## start from here to read in data
mean_sd_ldgr <- read.csv(paste(datpath, "mean_sd_ldgr.csv", sep = ""))

theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 16),
              axis.text = element_text(size = 16))

mean_sd_ldgr$N <- factor(mean_sd_ldgr$N, levels = c("lo","int","hi"))
mean_sd_ldgr$water <- factor(mean_sd_ldgr$water, levels = c("lo","hi"))
mean_sd_ldgr$invader_resident <- factor(mean_sd_ldgr$invader_resident, levels = c("pler_lapl","lapl_pler","pler_brho","brho_pler","lapl_brho","brho_lapl"))


levels(mean_sd_ldgr$invader_resident)= c("pler_lapl"=expression("*Plantago* invading *Layia*"),
                      "lapl_pler"=expression("*Layia* invading *Plantago*"),
                      "pler_brho"=expression("*Plantago* invading *Bromus*"),
                      "brho_pler"=expression("*Bromus* invading *Plantago*"),
                      "lapl_brho"=expression("*Layia* invading *Bromus*"),
                      "brho_lapl"=expression("*Bromus* invading *Layia*"))
                      
                      


ggplot(mean_sd_ldgr, aes(x=N, y=mean,ymin=mean-stdv,
                              ymax=mean+stdv, fill = water)) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  facet_wrap(~invader_resident,nrow=3,labeller = labeller(type=label_parsed)) +
  theme(strip.text.x = element_markdown())+
  ylab("Growth rate when rare") + xlab("N treatments") +
  scale_x_discrete(labels = c("Low","Intermediate","High")) +
  scale_fill_manual(name="Water treatments", labels = c("Dry","Wet"), 
                    values=c("azure3","azure4")) +
  geom_hline(yintercept = 0)+
  geom_errorbar(position = position_dodge(0.9),width=0.1)




