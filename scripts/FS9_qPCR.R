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
  ylab('log10 relative abundance of tet32 gene') +
  theme_bw()
tet32fig

tetWfig <- tetw %>% 
  ggplot(aes(x=Treatment, y=log10tetW, color=Day)) +
  scale_color_manual(values = c('7'='#E69F00', '11'='#CC0066', '14'='#999999')) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20)) +
  theme(axis.text.x=element_text(color = 'black', size = 14),
        axis.text.y=element_text(color = 'black', size=14)) + 
  ylab('log10 relative abundance of tetW gene') +
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
#Purpose: This code calculates AULC of days 0, 4, 7 for gene abundances of tetW, tet32, and aph2 for 
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

### Run one-way ANOVA
# One-way ANOVA test and multiple pairwise-comparisons test in R reference: 
# http://www.sthda.com/english/wiki/one-way-anova-test-in-r

# one-way ANOVA for tetW
colnames(sum_tetw)
tetw_anova <- aov(AULC~treatment, data=sum_tetw)
anova(tetw_anova)
# Results:
# Response: AULC
#           Df  Sum Sq  Mean Sq F value   Pr(>F)   
# treatment  2   14.878  7.4389  5.5632    0.009745 **
# Residuals 26   34.766  1.3372                    
# There is a significant difference in AULC between the treatment groups (group means are different). Which pairs of groups are different?

# Tukey multiple pairwise-comparisons
TukeyHSD(tetw_anova)
# Results:
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# Fit: aov(formula = AULC ~ treatment, data = sum_tetw)
# $treatment
#                     diff          lwr         upr         p adj
# INFinject-INFfeed   -1.51570154   -2.8359527  -0.1954504  0.0220943
# INFnm-INFfeed       -0.01680115   -1.3370523  1.3034500   0.9994489
# INFnm-INFinject     1.49890039    0.2138623   2.7839385   0.0198859

# Key for diff, lwr, upr, p adj:
# diff: difference between means of the two groups
# lwr, upr: the lower and the upper end point of the confidence interval at 95% (default)
# p adj: p-value after adjustment for the multiple comparisons.

# Pairwise t-test with corrections for multiple testing
pairwise.t.test(sum_tetw$AULC, sum_tetw$treatment, p.adjust.method = "BH")
# Pairwise comparisons using t tests with pooled SD 
# data:  sum_tetw$AULC and sum_tetw$treatment 
#           INFfeed INFinject
# INFinject 0.013   -        
# INFnm     0.975   0.013 

# Conclusion for tetW:
# Both Tukey's HSD and pairwise t-test agree that INFinject v INFfeed, and INFnm v INFinject have significantly different AULC

# one-way ANOVA for tet32
colnames(sum_tet32)
tet32_anova <- aov(AULC~treatment, data=sum_tet32)
anova(tet32_anova)
# Results:
# Response: AULC
#           Df  Sum Sq  Mean Sq F value   Pr(>F)   
# treatment  2  7.2275  3.6137  3.0171 0.06633 .
# Residuals 26  31.1421 1.1978                    
# The difference in AULC between the treatment groups approaches significant difference. Which pairs of groups are approaching significant difference?

# Tukey multiple pairwise-comparisons
TukeyHSD(tet32_anova)
# Results:
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# Fit: aov(formula = AULC ~ treatment, data = sum_tetw)
# $treatment
#                     diff        lwr         upr       p adj
# INFinject-INFfeed -1.09687068   -2.3464129  0.1526715 0.0933921
# INFnm-INFfeed     -0.09429518   -1.3438374  1.1552470 0.9808157
# INFnm-INFinject    1.00257550   -0.2136396  2.2187906 0.1207593

# Pairwise t-test with corrections for multiple testing
pairwise.t.test(sum_tet32$AULC, sum_tet32$treatment, p.adjust.method = "BH")
# Pairwise comparisons using t tests with pooled SD 
# data:  sum_tet32$AULC and sum_tet32$treatment 
#           INFfeed INFinject
# INFinject 0.076   -        
# INFnm     0.853   0.076 

# Conclusion for tet32:
# Both Tukey's HSD and pairwise t-test agree that there were no significant differences in AULC for tet32 gene abundance
# among the three treatment groups

# Graph b&w plot of AULC
#tetW
fig_tetw_aulc <- sum_tetw %>%ggplot(aes(x=treatment, y=AULC, color=treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20))+
  labs(y= 'AULC of relative tetW gene abundance', x= NULL) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme_bw()
fig_tetw_aulc

#tet32
fig_tet32_aulc <- sum_tet32 %>%ggplot(aes(x=treatment, y=AULC, color=treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20))+
  labs(y= 'AULC of relative tet32 gene abundance', x= NULL) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme_bw()
fig_tet32_aulc

#both tetW and tet32
fig_aulc <- plot_grid(fig_tet32_aulc, fig_tetw_aulc, labels = c('A', 'B'), label_size = 12)
fig_aulc
ggsave(fig_aulc,
       filename = './Fig8b_tet32tetW_AULC_INFfeedINFinjectINFnm.jpeg',
       width = 180,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

### Jules ###
sum_sal <- sal_data %>%
  arrange(time_point)  %>%
  group_by(pignum) %>%
  summarise(AULC=trapz(time_point, log_sal),
            treatment=unique(treatment))
