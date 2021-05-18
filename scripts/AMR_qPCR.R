#####################################################################################################
#FS9 tetW, tet32, aph2 qPCR - log10-fold gene abundance of INFinject and INFfeed relative to INFnm
#By Mou, KT; Stephens, A

#Purpose: This code graphs log10-fold gene abundances of tetW, tet32, and aph2 for INFinject, INFfeed, INFnm
#from qPCR studies

#Load library packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)

sessionInfo()
#R version 4.0.2 (2020-06-22)

########################################################################################################

tet32 <- read.csv('./data/FS9_tet32_qPCR_results.csv', stringsAsFactors = FALSE)
tetw <- read.csv('./data/FS9_tetW_qPCR_results.csv', stringsAsFactors = FALSE)
aph2 <- read.csv('./data/FS9_aph2_qPCR_results.csv', stringsAsFactors = FALSE)
tet32$Day <- factor(tet32$Day, levels=c("7", "11", "14"))
tetw$Day <- factor(tetw$Day, levels=c("7", "11", "14"))
aph2$Day <- factor(aph2$Day, levels=c("7", "11", "14"))
levels(tetw$Day) #"7"  "11" "14"
levels(tet32$Day) #"7"  "11" "14"
levels(aph2$Day) #"7"  "11" "14"

#Summary stats
tet32stats <- tet32 %>% group_by(Day,Treatment) %>% summarise(mean=mean(log10tet32), 
                                                                   n=n(), 
                                                                   sd=sd(log10tet32), 
                                                                   se=sd/n) #%>% write_csv('./results/tet32_qPCR_mean.csv')
tetWstats <- tetw %>% group_by(Day,Treatment) %>% summarise(mean=mean(log10tetW),
                                                            n=n(),
                                                            sd=sd(log10tetW),
                                                            se=sd/n) #%>% write_csv('./results/tetW_qPCR_mean.csv')
aph2stats <- aph2 %>% group_by(Day,Treatment) %>% summarise(mean=mean(log10aph2), 
                                                              n=n(), 
                                                              sd=sd(log10aph2), 
                                                              se=sd/n) #%>% write_csv('./results/aph2_qPCR_mean.csv')

########################################################################################################
#Purpose: This code calculates AULC of days 7, 11, 14 for gene abundances of tetW, tet32, and aph2 for 
#INFinject, INFfeed, INFnm from qPCR studies

library(tidyverse)
library(pracma)

#tet32
tet32$Day <- as.numeric(tet32$Day)
sapply(tet32,class) #make sure Day is numeric
sum_tet32 <- tet32 %>%
  arrange(Day)  %>%
  group_by(Pig) %>%
  summarise(AULC=trapz(Day, log10tet32),
            treatment=unique(Treatment))
write.csv(sum_tet32, file = "aulc_tet32.csv", col.names = TRUE)

#tetW
tetw$Day <- as.numeric(tetw$Day)
sapply(tetw,class) #make sure Day is numeric
sum_tetw <- tetw %>%
  arrange(Day)  %>%
  group_by(Pig) %>%
  summarise(AULC=trapz(Day, log10tetW),
            treatment=unique(Treatment))
write.csv(sum_tetw, file = "aulc_tetw.csv", col.names = TRUE)

#aph2
aph2$Day <- as.numeric(aph2$Day)
sapply(aph2,class) #make sure Day is numeric
sum_aph2 <- aph2 %>%
  arrange(Day)  %>%
  group_by(Pig.ID) %>%
  summarise(AULC=trapz(Day, log10aph2),
            treatment=unique(Treatment))
write.csv(sum_aph2, file = "aulc_aph2.csv", col.names = TRUE)

# Graph b&w plot of AULC
#tetW
(fig_tetw_aulc <- sum_tetw %>%ggplot(aes(x=treatment, y=AULC, color=treatment)) +
    scale_color_manual(values = c(INFinject='#00BA38',
                                  INFnm='#F8766D',
                                  INFfeed='#619CFF'),
                       name="Treatment") +
    geom_boxplot() + 
    geom_jitter(position=position_jitterdodge(jitter.width = .20))+
    labs(y= "AULC of tetW gene's relative abundance", x= NULL) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.position = "none") +
    ylim(5,10))

#tet32
(fig_tet32_aulc <- sum_tet32 %>%ggplot(aes(x=treatment, y=AULC, color=treatment)) +
    scale_color_manual(values = c(INFinject='#00BA38',
                                  INFnm='#F8766D',
                                  INFfeed='#619CFF')) +
    geom_boxplot() + 
    geom_jitter(position=position_jitterdodge(jitter.width = .20))+
    labs(y= "AULC of tet32 gene's relative abundance", x= NULL) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.position = "none") +
    ylim(5,10))

#aph2
(fig_aph2_aulc <- sum_aph2 %>%ggplot(aes(x=treatment, y=AULC, color=treatment)) +
    scale_color_manual(values = c(INFinject='#00BA38',
                                  INFnm='#F8766D',
                                  INFfeed='#619CFF'),
                       name="Treatment") +
    geom_boxplot() + 
    geom_jitter(position=position_jitterdodge(jitter.width = .20))+
    labs(y= "AULC of aph2 gene's relative abundance", x= NULL) +
    theme_bw() +
    theme(axis.text.y = element_text(size=12),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    ylim(5,10))

#all 3 figures together
(fig_aulc <- plot_grid(fig_tet32_aulc, fig_tetw_aulc, fig_aph2_aulc, labels = c('A', 'B', 'C'), label_size = 12, nrow=1, rel_widths = c(0.9, 0.9, 1.5)))
ggsave(fig_aulc,
       filename = './Fig6_tet32tetWaph2_AULC_INFfeedINFinjectINFnm.jpeg',
       width = 180,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

########################################################################################################
#Purpose: Create boxplots of relative log10-fold gene abundances of tetW, tet32, and aph2 for INFinject, INFfeed, INFnm

#Box plot figures (x-axis = treatment group)
tet32fig <- tet32 %>%
  ggplot(aes(x=Treatment, y=log10tet32, color=as.factor(Day))) +
  scale_color_manual(values = c('7'='#E69F00', '11'='#CC0066', '14'='#999999')) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
  theme(axis.text.x=element_text(color = 'black', size = 14),
        axis.text.y=element_text(color = 'black', size=14)) +
  ylab('log10 relative abundance of tet32 gene') +
  labs(color='Day')+
  theme_bw()
tet32fig+ theme(axis.title.x = element_blank())

tetWfig <- tetw %>%
  ggplot(aes(x=Treatment, y=log10tetW, color=as.factor(Day))) +
    scale_color_manual(values = c('7'='#E69F00', '11'='#CC0066', '14'='#999999')) +
    geom_boxplot() +
    geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
    theme(axis.text.x=element_text(color = 'black', size = 14),
          axis.text.y=element_text(color = 'black', size=4)) +
    ylab('log10 relative abundance of tetW gene') +
  labs(color='Day')+
    theme_bw()
tetWfig+ theme(axis.title.x = element_blank())

aph2fig <- aph2 %>%
  ggplot(aes(x=Treatment, y=log10aph2, color=as.factor(Day))) +
  scale_color_manual(values = c('7'='#E69F00', '11'='#CC0066', '14'='#999999')) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
  theme(axis.text.x=element_text(color = 'black', size = 30),
        axis.text.y=element_text(color = 'black', size=14))+
  ylab('log10 relative abundance of aph2 gene')+
  theme(axis.title.x = element_blank())+
  labs(color='Day')+
  theme_bw()
aph2fig+ theme(axis.title.x = element_blank())

#Box plot figures (x-axis = day)
tet32fig1 <- tet32 %>%
  ggplot(aes(x=Day, y=log10tet32, color=as.factor(Treatment))) +
  scale_color_manual(values = c('INFinject'='#00BA38', 'INFnm'='#F8766D', 'INFfeed'='#619CFF')) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
  theme(axis.text.x=element_text(color = 'black', size = 14),
        axis.text.y=element_text(color = 'black', size=14)) +
  ylab('log10 relative abundance of tet32 gene') +
  labs(color='Treatment')+
  theme_bw()
tet32fig1

tetWfig1 <- tetw %>%
  ggplot(aes(x=Day, y=log10tetW, color=as.factor(Treatment))) +
  scale_color_manual(values = c('INFinject'='#00BA38', 'INFnm'='#F8766D', 'INFfeed'='#619CFF')) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
  theme(axis.text.x=element_text(color = 'black', size = 14),
        axis.text.y=element_text(color = 'black', size=14)) +
  ylab('log10 relative abundance of tetW gene') +
  labs(color='Treatment')+
  theme_bw()
tetWfig1

aph2fig1 <- aph2 %>%
  ggplot(aes(x=Day, y=log10aph2, color=as.factor(Treatment))) +
  scale_color_manual(values = c('INFinject'='#00BA38', 'INFnm'='#F8766D', 'INFfeed'='#619CFF')) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
  theme(axis.text.x=element_text(color = 'black', size = 14),
        axis.text.y=element_text(color = 'black', size=14)) +
  ylab('log10 relative abundance of aph2 gene') +
  labs(color='Treatment')+
  
  theme_bw()
aph2fig1

fig5 <- plot_grid(tet32fig+ theme(axis.title.x = element_blank()), tetWfig+ theme(axis.title.x = element_blank()),aph2fig+ theme(axis.title.x = element_blank()), tet32fig1, tetWfig1, aph2fig1, labels = c('A','B','C','D','E','F'), label_size = 12)
fig5
ggsave(fig5,
       filename = './Fig5_tet32tetWaph2_1_qPCR_INFfeedINFinjectINFnm.jpeg',
       width = 250,
       height = 180,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')