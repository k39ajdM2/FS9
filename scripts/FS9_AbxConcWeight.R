#####################################################################################################
#FS9 Antibiotic concentration and weight
#Kathy Mou

#Purpose: This code graphs concentration of oxytetracycline of each tissue for each group and also relative to weight.
#This was taken from Jules' code (FS1_ABX_conc.R) with some modifications to fit FS9 dataset.

#Load library packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)

sessionInfo()
#R version 3.6.3 (2020-02-29)

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
