#####################################################################################################
#FS9 Antibiotic concentration, weight, average daily gain (ADG), lung lesion
#By Mou, KT

#Purpose: This code graphs concentration of oxytetracycline in each tissue for each group and also relative to weight, ADG, and lung lesion
#This was taken from Jules' code (FS1_ABX_conc.R) with some modifications to fit FS9 dataset.

#Load library packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
#library("ggsignif")

sessionInfo()
#R version 4.0.2 (2020-06-22)

########################################################################################################

#Import file
tis <- read.csv('./data/FS9_WeightOxytet.csv', stringsAsFactors = FALSE)  # reads in data, already cleaned a little

#Must organize "tis" file in order to run the various scripts below
tis$Day <- as.numeric(gsub('D', '', tis$Day))             # replaces the 'D' in the time column with '' (nothing)
tis.melt <- melt(tis, id.vars = c(1:4), measure.vars = c(5:7))          # converts to long dataframe format for easy plotting
tis.melt                                                                # just checking on the new long dataframe
tis.melt$value <- gsub('NF', 0, tis.melt$value)         # replaces 'NF' with 0
tis.melt$value <- as.numeric(tis.melt$value)            # forces the value column to be numeric 
tis.melt <- na.exclude(tis.melt)                                        # removes NAs which were introduced in the previous line
tis.melt$DayXTreatment <- paste(tis.melt$Day, tis.melt$Treatment, sep = '_') # creates a 'DayXTreatment column', it's just the 'Day' and 'Treatment' columns pasted together
tis.melt$Day <- factor(tis.melt$Day)                                     # makes the 'Day' column a factor (this means it is categorical data, not continuous data)
tis.melt$Tissue <- ifelse(tis.melt$variable == 'PlasmaOxytet', 'Plasma', 
                          ifelse(tis.melt$variable == 'LungOyxtet', 'Lung',
                                 ifelse(tis.melt$variable == 'NasalOxytet', 'Nasal', NA)))   #Replace original names in 'variable' column with with new names
unique(tis.melt$Treatment)
tis.melt$Treatment <- factor(tis.melt$Treatment, levels = c('NONINFnm', 'INFnm', 'INFinject', 'INFfeed'))

#######################################################################################################################################
#Oxytetracycline Levels, Weight (Figure 3)

#Plasma oxytetracycline only
plasmaoxytet <- tis.melt %>% filter(Tissue == 'Plasma') %>% filter(Treatment %in% c('INFinject', "INFfeed")) %>% 
  ggplot(aes(x=Day, y=value, group=DayXTreatment, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                               INFfeed='#619CFF')) +
  geom_boxplot() + 
  ylab('Concentration of Oxytet \n (ng/mL plasma)') + 
  theme_bw() + 
  geom_jitter(position=position_jitterdodge(jitter.width = 1)) +
  ylim(0,800)
plasmaoxytet

#Lung oxytetracycline only
lungoxytet <- tis.melt %>% filter(Tissue == 'Lung') %>% filter(Treatment %in% c('INFinject', "INFfeed")) %>% 
  ggplot(aes(x=Day, y=value, group=DayXTreatment, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                               INFfeed='#619CFF')) +
  geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 1)) +
  ylab('Concentration of Oxytet \n (ng/2g lung tissue)') + 
  theme_bw() +
  ylim(0,450)
lungoxytet

#Nasal oxytetracycline only
nasaloxytet <- tis.melt %>% filter(Tissue == 'Nasal') %>% filter(Treatment %in% c('INFinject', 'INFfeed')) %>% 
  ggplot(aes(x=Day, y=value, group=DayXTreatment, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                               INFfeed='#619CFF')) +
  geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 1)) +
  ylab('Concentration of Oxytet \n (ng/mL nasal wash)') +
  theme_bw() +
  ylim(0,550)
nasaloxytet

#Oxytet stats
tis.melt$group <- paste(tis.melt$Pig, tis.melt$Day, tis.melt$Tissue, sep = '_')
Lung.abx <- tis.melt %>% filter(Tissue == 'Lung') #subset lung data
tmp <- tis.melt %>% group_by(Day, Tissue, Treatment) %>% summarise(mean=mean(value), 
                                                                   n=n(), 
                                                                   sd=sd(value), 
                                                                   se=sd/n) %>% write_csv('./results/mean_abx_conc.csv')

#Weight and oxytetracycline plots

#Convert weight from character to numeric in order to run linear regression with Weight as x-axis
class(tis.melt$Weight) #character
tis.melt$Weight <- as.numeric(tis.melt$Weight)
class(tis.melt$Weight) #numeric

#Concentration of oxytetracycline relative to weight, plasma, days 11 and 14
plasmaweight <- tis.melt %>% 
  filter(Tissue=="Plasma") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF',
                                NONINFnm="#C77CFF")) +
  geom_smooth(method = 'lm') + 
  geom_point() + 
  theme_bw() +
  ylab('Concentration of Oxytet \n (ng/mL plasma)') + 
  xlab("Weight (lbs)") +
  facet_wrap(vars(Day), scales = "free") +
  xlim(13,34) +
  ylim(0,800)
plasmaweight

#Concentration of oxytetracycline relative to weight, lung, days 11 and 14
lungweight <- tis.melt %>% 
  filter(Tissue=="Lung") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF',
                                NONINFnm="#C77CFF")) +
  geom_smooth(method = 'lm') + 
  geom_point() + 
  theme_bw() +
  ylab('Concentration of Oxytet \n (ng/2g lung tissue)') + 
  xlab("Weight (lbs)") +
  facet_wrap(vars(Day), scales = "free") +
  xlim(13,34) +
  ylim(0,450)
lungweight

#Concentration of oxytetracycline relative to weight, nasal, days 11 and 14
nasalweight <- tis.melt %>% 
  filter(Tissue=="Nasal") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF',
                                NONINFnm="#C77CFF")) +
  geom_smooth(method = 'lm') + 
  geom_point() + 
  theme_bw() +
  ylab('Concentration of Oxytet \n (ng/mL nasal wash)') + 
  xlab("Weight (lbs)") +
  theme(axis.title.y = element_text(size=10)) +
  facet_wrap(vars(Day), scales = "free") +
  xlim(13,34) +
  ylim(0,550)
nasalweight

#Combine figures
fig_3 <- plot_grid(lungoxytet, lungweight, nasaloxytet, nasalweight, plasmaoxytet, plasmaweight, 
                   labels = c('A', 'B', "C", "D", "E", "F"), 
                   label_size = 12, 
                   ncol = 2,
                   vjust = 1.4,
                   rel_widths = c(1, 1.5))
fig_3

ggsave("Fig3_OxytetLevelsWeight.tiff", plot=fig_3, width = 10, height = 8, dpi = 200, units =c("in"))



########################################################################################################
#ADG (Figure 2)

#Import file
adg <- read.csv('./data/FS9_AverageDailyGain.csv', stringsAsFactors = FALSE)
colnames(adg)
adg1 <- pivot_longer(adg, cols=c("D0", "D11", "D14"), names_to="Day", values_to="Weight_lbs")
adg2 <- pivot_longer(adg1, cols=c("D11_ADG", "D14_ADG"), names_to="Day_ADG", values_to="ADG")
adg2$Day_ADG <- as.character((adg2$Day_ADG))
adg2$Day <- as.character((adg2$Day))

stats <- adg2 %>% 
  group_by(Treatment, Day_ADG) %>% 
  summarise(Sum=sum(ADG, na.rm = TRUE), Mean=(mean(ADG, na.rm = TRUE)), sd = sd(ADG, na.rm = TRUE)) %>% 
write.csv(stats, file= "FS9_ADGstats.csv")

stats2 <- adg %>% 
  group_by(Treatment) %>% 
  summarise(Sum=sum(D0, na.rm = TRUE), Mean=(mean(D0, na.rm = TRUE)), sd = sd(D0, na.rm = TRUE))
write.csv(stats2, file="FS9_D0weight.csv")

colnames(adg)

stats3 <- adg %>% 
  group_by(Treatment) %>% 
  mutate(two_week_wt = (D14 - D0)) %>% 
  summarise(Sum=sum(two_week_wt, na.rm = TRUE), Mean=(mean(two_week_wt, na.rm= TRUE)), sd = sd(two_week_wt, na.rm = TRUE))
write.csv(stats3, file="FS9_two_week_wt.csv")

#Figure 2a
fig_day0adg <- adg %>% 
  ggplot(aes(x=Treatment, y=D0, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF',
                                NONINFnm="#C77CFF")) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20))+
  labs(y= 'Weight (lbs) on day 0', x= NULL) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme_bw() + 
  theme(legend.position = "none")
fig_day0adg

#Figure 2b
fig_day11adg <- adg %>% 
  ggplot(aes(x=Treatment, y=D11_ADG, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF',
                                NONINFnm="#C77CFF")) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20))+
  labs(y= 'Average daily gain (lbs) from day 0 to day 11', x= NULL) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(0,1.2)
fig_day11adg

#Figure 2c
fig_day14adg <- adg %>% 
  ggplot(aes(x=Treatment, y=D14_ADG, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF',
                                NONINFnm="#C77CFF")) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20))+
  labs(y= 'Average daily gain (lbs) from day 0 to day 14', x= NULL) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme_bw() + 
  ylim(0,1.2)
fig_day14adg

#Combined Weight + ADG figures
fig_adg3plots <- plot_grid(fig_day0adg, fig_day11adg, fig_day14adg, labels = c('A', 'B', 'C'), label_size = 12, nrow = 1, rel_widths = c(1, 1, 1.2))
fig_adg3plots
ggsave("Fig2_ADG_abc.tiff", plot=fig_adg3plots, width = 15, height = 5, dpi = 500, units =c("in"))

########################################################################################################
#Lung Lesion (Supplemental Figure)

lung <- read.csv('./data/FS9_LungLesionScores.csv', stringsAsFactors = FALSE)
colnames(lung)
fig_lung <- lung %>% 
  ggplot(aes(x=Day, y=WeightedAverage, color=Treatment)) +
  scale_color_manual(values = c(NONINFnm="#C77CFF",
                                INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + 
  geom_jitter(aes(color= Treatment), position=position_jitterdodge(jitter.width = .20))+
  labs(y= 'Lung Lesion Weighted Average Score', x= NULL) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme_bw()
fig_lung

stats <- lung %>% 
  group_by(lung$Treatment, Day) %>% 
  summarise(Sum=sum(WeightedAverage), Mean=(mean(WeightedAverage)), sd = sd(WeightedAverage), se = sd(WeightedAverage)/sqrt(n()))
write.csv(stats, file= "FS9_LungLesionStats.csv")

ggsave("SuppFig1_LungLesionSeverity.tiff", plot=fig_lung, width = 5, height = 7, dpi = 500, units =c("in"))
