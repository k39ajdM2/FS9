#####################################################################################################
#FS9 tetW, tet32, aph2 qPCR - log10-fold gene abundance of INFinject and INFfeed relative to INFnm
#Kathy Mou

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
tet32$Day <- factor(tet32$Day, levels=c("7", "11", "14"))
tetw$Day <- factor(tetw$Day, levels=c("7", "11", "14"))
levels(tetw$Day) #"7"  "11" "14"
levels(tet32$Day) #"7"  "11" "14"

str(tet32) #see 

#Two-way ANOVA with repeated measures, Tukey multiple pairwise-comparisons: 
#http://www.sthda.com/english/wiki/two-way-anova-test-in-r#multiple-pairwise-comparison-between-the-means-of-groups

#Two way ANOVA
tet32aov <- aov(log10tet32 ~ Treatment * Day, data=tet32) 
summary(tet32aov)
#               Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment      2   3.42   1.711   2.383    0.0989 .  
#Day            2  20.19  10.093  14.058    6.11e-06 ***
#Treatment:Day  4   8.41   2.103   2.930    0.0260 *  
#Residuals     78  56.00   0.718                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#These results would lead us to believe that the levels of Day, and Day + Treatment interaction
#are associated with significant different mean log10tet32 abundances.
#In addition, the interaction between Day and Treatment indicates that the relationship between mean log10tet32 abundances
#and Treatment depend on the day sampled.

tetWaov <- aov(log10tetW ~ Treatment * Day, data=tetw)
summary(tetWaov)
#```````````````Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment      2   6.88   3.440   5.110    0.00821 ** 
#Day            2  18.05   9.025  13.408    9.88e-06 ***
#Treatment:Day  4  11.16   2.791   4.146    0.00425 ** 
#Residuals     78  52.50   0.673                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#These results would lead us to believe that the levels of Day or Treatment, and Day + Treatment interaction
#are associated with significant different mean log10tet32 abundances.
#In addition, the interaction between Day and Treatment indicates that the relationship between mean log10tet32 abundances
#and Day depend on the type of Treatment administered.

#Tukey multiple pairwise-comparisons
TukeyHSD(tet32aov, which = "Treatment:Day")
#can't save as csv so I copied the results in console to 
#tet32_qPCR_summarystats_Tukey.xlsx
TukeyHSD(tetWaov, which = "Treatment:Day")
#can't save as csv so I copied the results in console to
#tetW_qPCR_summarystats_Tukey.xlsx

#Summary stats
tet32stats <- tet32 %>% group_by(Day,Treatment) %>% summarise(mean=mean(log10tet32), 
                                                                   n=n(), 
                                                                   sd=sd(log10tet32), 
                                                                   se=sd/n) %>% write_csv('./results/tet32_qPCR_mean.csv')

tetWstats <- tetw %>% group_by(Day,Treatment) %>% summarise(mean=mean(log10tetW),
                                                            n=n(),
                                                            sd=sd(log10tetW),
                                                            se=sd/n) %>% write_csv('./results/tetW_qPCR_mean.csv')

tet32fig <- tet32 %>% 
  ggplot(aes(x=Treatment, y=log10tet32, color=Day)) +
  scale_color_manual(values = c('7'='#E69F00', '11'='#CC0066', '14'='#999999')) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
  theme(axis.text.x=element_text(color = 'black', size = 14),
        axis.text.y=element_text(color = 'black', size=14)) + 
  ylab('log10-fold change of tet32 gene abundance') +
  theme_bw()
tet32fig

tetWfig <- tetw %>% 
  ggplot(aes(x=Treatment, y=log10tetW, color=Day)) +
  scale_color_manual(values = c('7'='#E69F00', '11'='#CC0066', '14'='#999999')) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
  theme(axis.text.x=element_text(color = 'black', size = 14),
        axis.text.y=element_text(color = 'black', size=14)) + 
  ylab('log10-fold change of tetW gene abundance') +
  theme_bw()
tetWfig

#Combine tet32 and tetW graphs and save combined figure
fig8 <- plot_grid(tet32fig, tetWfig, labels = c('A', 'B'), label_size = 12)
fig8
ggsave(fig8,
       filename = './Fig8_tet32tetW_qPCR_INFfeedINFinjectINFnm.jpeg',
       width = 180,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

########################################################################################################
#Purpose: This code calculates AUC of days 0, 4, 7 for gene abundances of tetW, tet32, and aph2 for INFinject, INFfeed, INFnm
#from qPCR studies

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

#tetW
tetw$Day <- as.numeric(tetw$Day)
sapply(tetw,class) #make sure Day is numeric
sum_tetw <- tetw %>%
  arrange(Day)  %>%
  group_by(Pig) %>%
  summarise(AULC=trapz(Day, log10tetW),
            treatment=unique(Treatment))

## graph AULC data, run one-way ANOVA


### Jules ###
sum_sal <- sal_data %>%
  arrange(time_point)  %>%
  group_by(pignum) %>%
  summarise(AULC=trapz(time_point, log_sal),
            treatment=unique(treatment))
