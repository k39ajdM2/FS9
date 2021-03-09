#####################################################################################################
#FS9 Antibiotic concentration, weight, average daily gain, Bordetella bronchiseptica and Pasteurella multocida colonization, lung lesion
#Kathy Mou

#Purpose: This code graphs concentration of oxytetracycline of each tissue for each group and also relative to weight; ADG, B. bronchiseptica and
#P. multocida colonization in each tissue; lung lesion
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

######################################################################################################################################################
#Oxytetracycline Levels, Weight

#Oxytetracycline concentration in INFinject and INFfeed in all tissue samples
fig2 <- tis.melt %>% filter(Treatment %in% c('INFinject', 'INFfeed')) %>%
  ggplot(aes(x=Day, y=value, group=DayXTreatment, color=Treatment)) +
  geom_boxplot() +
  ylab('Concentration of Oxytetracycline (ng/mL)') + xlab('Day') + scale_y_log10(labels=scales::scientific) +
  facet_wrap(~Tissue)+
  geom_jitter(position=position_jitterdodge(jitter.width = .20))+
  theme_bw() + scale_color_manual(values = c(INFinject='#00BA38', INFfeed='#619CFF'))
fig2
ggsave(fig2,
       filename = './Fig_OxytetLevels_LungNasalPlasma_INFfeedINFinject.jpeg',
       width = 180,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

class(tis.melt$Day)
#tis.melt$Day <- as.numeric(as.character(tis.melt$Day))
#tis.melt$Weight <- as.numeric(as.character(tis.melt$Weight))

#Plasma oxytetracycline only
plasmaoxytet <- tis.melt %>% filter(Tissue == 'Plasma') %>% filter(Treatment %in% c('INFinject', "INFfeed")) %>% 
  ggplot(aes(x=Day, y=value, group=DayXTreatment, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                               INFfeed='#619CFF')) +
  geom_boxplot() + 
  ylab('Concentration of Oxytet \n (ng/mL plasma)') + 
  theme_bw() + 
  geom_jitter(position=position_jitterdodge(jitter.width = 1))
plasmaoxytet

#Lung oxytetracycline only
lungoxytet <- tis.melt %>% filter(Tissue == 'Lung') %>% filter(Treatment %in% c('INFinject', "INFfeed")) %>% 
  ggplot(aes(x=Day, y=value, group=DayXTreatment, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                               INFfeed='#619CFF')) +
  geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 1)) +
  ylab('Concentration of Oxytet \n (ng/2g lung tissue)') + 
  theme_bw()
lungoxytet

#Nasal oxytetracycline only
nasaloxytet <- tis.melt %>% filter(Tissue == 'Nasal') %>% filter(Treatment %in% c('INFinject', 'INFfeed')) %>% 
  ggplot(aes(x=Day, y=value, group=DayXTreatment, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                               INFfeed='#619CFF')) +
  geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 1)) +
  ylab('Concentration of Oxytet \n (ng/mL nasal wash)') +
  theme_bw()
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
  xlim(13,34)
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
  xlim(13,34)
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
  xlim(13,34)
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
#Colonization

# install.packages("mdthemes")
library("mdthemes")

#Import file
bb <- read.csv('./data/FS9_Bordetella_bronchiseptica_colonization.csv', stringsAsFactors = FALSE)  # reads in data, already cleaned a little
pm <- read.csv('./data/FS9_Pasteurella_multocida_colonization.csv', stringsAsFactors = FALSE)

#Bordetella
bb2 <- pivot_longer(bb, cols=c("Nasal._Wash", "Tonsil", "Tracheal_Wash", "Lung_Lavage", "Lung_Tissue"), names_to="Sample", values_to="log10CFU/g_or_mL")
bb2$Sample <- as.character((bb2$Sample))
bb2$Sample[bb2$Sample == "Nasal._Wash"] <- "Nasal_Wash"
colnames(bb2)
names(bb2)[names(bb2)== "log10CFU/g_or_mL"] <- "log10CFU"

#Samples
bb2 %>% filter(Sample == 'Nasal_Wash') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per mL in nasal wash')

bb2 %>% filter(Sample == 'Tonsil') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per gram in tonsil')

bb2 %>% filter(Sample == 'Tracheal_Wash') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per mL in tracheal wash')

bb2 %>% filter(Sample == 'Lung_Lavage') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per mL in lung lavage')

bb2 %>% filter(Sample == 'Lung_Tissue') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per gram in lung tissue')

stats <- bb2 %>% 
  group_by(Treatment, Day, Sample) %>% 
  summarise(Sum=sum(log10CFU, na.rm = TRUE), Mean=(mean(log10CFU, na.rm = TRUE)), sd = sd(log10CFU, na.rm = TRUE))
stats$sd <- round(stats$sd, 3)
stats$Mean <- round(stats$Mean, 3)
stats$meanSD <- with(stats, paste0(Mean, sep="+/-", sd)) 
write.csv(stats, file= "FS9_Q2_Bordetella_colonizationstats.csv")

#Pasteurella
pm2 <- pivot_longer(pm, cols=c("Nasal_Wash", "Tonsil"), names_to="Sample", values_to="log10CFU")
pm2$Sample <- as.character((pm2$Sample))
colnames(pm2)

#Nasal wash
pm_nasal <- pm2 %>% filter(Sample == 'Nasal_Wash') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + labs(y= '*Pasteurella multocida* log10CFU per mL in nasal wash') + md_theme_linedraw() +
  ylim(0, 7)
pm_nasal
ggsave(pm_nasal,
       filename = './figure_Pm_nasal.jpeg',
       width = 120,
       height = 140,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

#Tonsil
pm2 %>% filter(Sample == 'Tonsil') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_boxplot() + ylab('Pasteurella multocida log10CFU per gram in tonsil') + md_theme_linedraw()

stats2 <- pm2 %>% 
  group_by(Treatment, Day, Sample) %>% 
  summarise(Sum=sum(log10CFU, na.rm = TRUE), Mean=(mean(log10CFU, na.rm = TRUE)), sd = sd(log10CFU, na.rm = TRUE))
stats2$sd <- round(stats2$sd, 3)
stats2$Mean <- round(stats2$Mean, 3)
stats2$meanSD <- with(stats2, paste0(Mean, sep=" + ", sd)) 
write.csv(stats2, file= "FS9_Q2_Pasteurella_colonizationstats.csv")

########################################################################################################
#ADG

#Import file
adg <- read.csv('./data/FS9_AverageDailyGain.csv', stringsAsFactors = FALSE)
colnames(adg)
adg1 <- pivot_longer(adg, cols=c("D0", "D11", "D14"), names_to="Day", values_to="Weight_lbs")
adg2 <- pivot_longer(adg1, cols=c("D11_ADG", "D14_ADG"), names_to="Day_ADG", values_to="ADG")
adg2$Day_ADG <- as.character((adg2$Day_ADG))
adg2$Day <- as.character((adg2$Day))

stats <- adg2 %>% 
  group_by(Treatment, Day_ADG) %>% 
  summarise(Sum=sum(ADG, na.rm = TRUE), Mean=(mean(ADG, na.rm = TRUE)), sd = sd(ADG, na.rm = TRUE))
write.csv(stats, file= "FS9_ADGstats.csv")

#Old figure 1
fig_adg <- adg2 %>% 
  ggplot(aes(x=Day_ADG, y=ADG, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF',
                                NONINFnm="#C77CFF")) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20))+
  labs(y= 'Average Daily Gain (pounds)', x= NULL) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  scale_x_discrete(breaks=c("D4_ADG", "D7_ADG"),
                   labels=c("D4", "D7")) +
  theme_bw()
fig_adg

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
  theme(legend.position = "none")
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
  theme_bw()
fig_day14adg

#Combined Weight + ADG figures
fig_adg3plots <- plot_grid(fig_day0adg, fig_day11adg, fig_day14adg, labels = c('A', 'B', 'C'), label_size = 12, nrow = 1, rel_widths = c(1, 1, 1.2))
fig_adg3plots
ggsave("Fig2_ADG_abc.tiff", plot=fig_adg3plots, width = 15, height = 5, dpi = 500, units =c("in"))

#ggsave(fig_adg,
#        filename = './Fig2_ADG_All.jpeg',
#       width = 160,
#       height = 200,
#       device = 'jpeg',
#       dpi = 300,
#       units = 'mm')

########################################################################################################
#Lung Lesion

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

#Attempting to add significance bar with geom_signif package, but the line and asterisk 
#keeps coming out as blue color. Not sure how to make it a different color.
#fig_lung + geom_signif(
#  y_position = c(7.5), xmin = c(2.1), xmax = c(2.4),
#  annotation = c("*"), tip_length = 0, textsize = 5)

stats <- lung %>% 
  group_by(lung$Treatment, Day) %>% 
  summarise(Sum=sum(WeightedAverage), Mean=(mean(WeightedAverage)), sd = sd(WeightedAverage), se = sd(WeightedAverage)/sqrt(n()))
write.csv(stats, file= "FS9_LungLesionStats.csv")

ggsave("SuppFig1_LungLesionSeverity.tiff", plot=fig_lung, width = 5, height = 7, dpi = 500, units =c("in"))

#ggsave(fig_lung,
#       filename = './figure_lung.jpeg',
#       width = 160,
#       height = 200,
#       device = 'jpeg',
#       dpi = 300,
#       units = 'mm')
