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

sessionInfo()
#R version 4.0.2 (2020-06-22)

########################################################################################################

#Import file
tis <- read.csv('./data/FS9_WeightOxytet.csv', stringsAsFactors = FALSE)  # reads in data, already cleaned a little

#Organize "tis" file
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

#Oxytetracycline concentration in two groups in all tissue samples
fig2 <- tis.melt %>% filter(Treatment %in% c('INFinject', 'INFfeed')) %>%
  ggplot(aes(x=Day, y=value, group=DayXTreatment, fill=Treatment)) +
  geom_boxplot() +
  ylab('Concentration of Oxytetracycline (ng/mL)') + xlab('Day') + scale_y_log10(labels=scales::scientific) +
  facet_wrap(~Tissue, scales = 'free')+
  theme_bw() + scale_fill_manual(values = c('#E69F00', '#999999'))
fig2
ggsave(fig2,
       filename = './figure2.jpeg',
       width = 180,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

class(tis.melt$Day)
tis.melt$Day <- as.numeric(as.character(tis.melt$Day))
tis.melt$Weight <- as.numeric(as.character(tis.melt$Weight))

#Oxytet concentration in all four groups, all tissues
tis.melt %>% 
  ggplot(aes(x=Day, y=value, fill=Treatment, color=Treatment)) +
  geom_point(shape=21, alpha=.5) + geom_smooth() +
  ylab('concentration ng/mL') + scale_y_log10() +
  facet_wrap(~Tissue, scales = 'free')+
  theme_bw() + scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999'))+
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) 

#BW plot of oxytet concentration, all four groups, all tissues
tis.melt %>% 
  ggplot(aes(x=Day, y=value, group=DayXTreatment, fill=Treatment)) +
  geom_boxplot() +
  ylab('concentration ng/mL') +
  facet_wrap(~Tissue, scales = 'free')+
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  theme_bw()

#Plasma oxytetracycline only
tis.melt %>% filter(Tissue == 'Plasma') %>%
  ggplot(aes(x=Day, y=value, group=DayXTreatment, fill=Treatment)) +
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_boxplot() + ylab('concentration ng/mL') + ggtitle('Oxytet concentrations: Plasma')

#Lung oxytetracycline only
tis.melt %>% filter(Tissue == 'Lung') %>%
  ggplot(aes(x=Day, y=value, group=DayXTreatment, fill=Treatment)) +
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_boxplot() + ylab('concentration ng/mL') + ggtitle('Oxytet concentrations: Lung')

#Nasal oxytetracycline only
tis.melt %>% filter(Tissue == 'Nasal') %>%
  ggplot(aes(x=Day, y=value, group=DayXTreatment, fill=Treatment)) +
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_boxplot() + ylab('concentration ng/mL') + ggtitle('Oxytet concentrations: Nasal')


tis.melt$group <- paste(tis.melt$Pig, tis.melt$Day, tis.melt$Tissue, sep = '_')
Lung.abx <- tis.melt %>% filter(Tissue == 'Lung') #subset lung data
tmp <- tis.melt %>% group_by(Day, Tissue, Treatment) %>% summarise(mean=mean(value), 
                                                                   n=n(), 
                                                                   sd=sd(value), 
                                                                   se=sd/n) %>% write_csv('./results/mean_abx_conc.csv')

#Concentration of oxytetracycline relative to weight, plasma, day 4
fig1_d4 <- tis.melt %>% 
  filter(Day==4) %>%
  filter(Tissue=="Plasma") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_smooth(method = 'lm') + geom_point() + theme_bw() +ylab('Concentration of Oxytetracycline (ng/mL)') + xlab("Weight (kg)")
fig1_d4

#Concentration of oxytetracycline relative to weight, plasma, day 7
fig1_d7 <- tis.melt %>% 
  filter(Day==7) %>%
  filter(Tissue=="Plasma") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_smooth(method = 'lm') + geom_point() + theme_bw() +ylab('Concentration of Oxytetracycline (ng/mL)') + xlab("Weight (kg)")
fig1_d7

#Combine figures
fig1_plasma <- plot_grid(fig1_d4, fig1_d7, labels = c('A', 'B'), label_size = 12)

#Save combined figure
ggsave(fig1_plasma,
       filename = './results/figure1_plasma.jpeg',
       width = 180,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')



#Concentration of oxytetracycline relative to weight, lung, day 4
fig1_d4_lung <- tis.melt %>% 
  filter(Day==4) %>%
  filter(Tissue=="Lung") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_smooth(method = 'lm') + geom_point() + theme_bw() +ylab('Concentration of Oxytetracycline (ng/mL)') + xlab("Weight (kg)")
fig1_d4_lung

#Concentration of oxytetracycline relative to weight, lung, day 7
fig1_d7_lung <- tis.melt %>% 
  filter(Day==7) %>%
  filter(Tissue=="Lung") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_smooth(method = 'lm') + geom_point() + theme_bw() +ylab('Concentration of Oxytetracycline (ng/mL)') + xlab("Weight (kg)")
fig1_d7_lung

#Combine figures
fig1_lung <- plot_grid(fig1_d4_lung, fig1_d7_lung, labels = c('A', 'B'), label_size = 12)
fig1_lung

#Save combined figure
ggsave(fig1_lung,
       filename = './results/figure1_lung.jpeg',
       width = 180,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')


#Concentration of oxytetracycline relative to weight, nasal, day 4
fig1_d4_nasal <- tis.melt %>% 
  filter(Day==4) %>%
  filter(Tissue=="Nasal") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_smooth(method = 'lm') + geom_point() + theme_bw() +ylab('Concentration of Oxytetracycline (ng/mL)') + xlab("Weight (kg)")
fig1_d4_nasal

#Concentration of oxytetracycline relative to weight, nasal, day 7
fig1_d7_nasal <- tis.melt %>% 
  filter(Day==7) %>%
  filter(Tissue=="Nasal") %>% 
  ggplot(aes(x=Weight, y=value, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', NONINFnm='#56B4E9', INFfeed='#999999')) +
  geom_smooth(method = 'lm') + geom_point() + theme_bw() +ylab('Concentration of Oxytetracycline (ng/mL)') + xlab("Weight (kg)")
fig1_d7_nasal

#Combine figures
fig1_nasal <- plot_grid(fig1_d4_nasal, fig1_d7_nasal, labels = c('A', 'B'), label_size = 12)
fig1_nasal

#Save combined figure
ggsave(fig1_nasal,
       filename = './results/figure1_nasal.jpeg',
       width = 180,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

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
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', INFfeed='#999999')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per mL in nasal wash')

bb2 %>% filter(Sample == 'Tonsil') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', INFfeed='#999999')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per gram in tonsil')

bb2 %>% filter(Sample == 'Tracheal_Wash') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', INFfeed='#999999')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per mL in tracheal wash')

bb2 %>% filter(Sample == 'Lung_Lavage') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', INFfeed='#999999')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per mL in lung lavage')

bb2 %>% filter(Sample == 'Lung_Tissue') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', INFfeed='#999999')) +
  geom_boxplot() + ylab('Bordetella bronchiseptica log10CFU per gram in lung tissue')


#Pasteurella
pm2 <- pivot_longer(pm, cols=c("Nasal_Wash", "Tonsil"), names_to="Sample", values_to="log10CFU")
pm2$Sample <- as.character((pm2$Sample))
colnames(pm2)

#Nasal wash
pm_nasal <- pm2 %>% filter(Sample == 'Nasal_Wash') %>%
  ggplot(aes(x=Day, y=log10CFU, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', INFfeed='#999999')) +
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
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', INFfeed='#999999')) +
  geom_boxplot() + ylab('Pasteurella multocida log10CFU per gram in tonsil')


########################################################################################################
#ADG

#Import file
adg <- read.csv('./data/FS9_AverageDailyGain.csv', stringsAsFactors = FALSE)
colnames(adg)
adg1 <- pivot_longer(adg, cols=c("DNEG8", "D4", "D7"), names_to="Day", values_to="Weight_kg")
adg2 <- pivot_longer(adg1, cols=c("D4_ADG", "D7_ADG"), names_to="Day_ADG", values_to="ADG")
adg2$Day_ADG <- as.character((adg2$Day_ADG))
adg2$Day <- as.character((adg2$Day))

stats <- adg2 %>% 
  group_by(adg2$Treatment, Day_ADG) %>% 
  summarise(Sum=sum(ADG, na.rm = TRUE), Mean=(mean(ADG, na.rm = TRUE)), sd = sd(ADG, na.rm = TRUE))
write.csv(stats, file= "ADGstats.csv")

fig_adg <- adg2 %>% 
  ggplot(aes(x=Day_ADG, y=ADG, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00',
                                INFnm='#CC0066',
                                INFfeed='#999999',
                                NONINFnm="#56B4E9" )) +
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

ggsave("ADG.tiff", plot=fig_adg, width = 5, height = 7, dpi = 500, units =c("in"))

#ggsave(fig_adg,
#        filename = './figure_ADG.jpeg',
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
  scale_color_manual(values = c(INFinject='#E69F00',
                                INFnm='#CC0066',
                                INFfeed='#999999',
                                NONINFnm="#56B4E9" )) +
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(jitter.width = .20))+
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
  summarise(Sum=sum(WeightedAverage), Mean=(mean(WeightedAverage)), sd = sd(WeightedAverage))
write.csv(stats, file= "lunglesionstats.csv")

ggsave("LungLesionSeverity.tiff", plot=fig_lung, width = 5, height = 7, dpi = 500, units =c("in"))

#ggsave(fig_lung,
#       filename = './figure_lung.jpeg',
#       width = 160,
#       height = 200,
#       device = 'jpeg',
#       dpi = 300,
#       units = 'mm')
