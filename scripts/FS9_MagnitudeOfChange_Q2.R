#####################################################################################################
#FS9 Magnitude of Change - Infected: nm vs feed or inject
#Kathy Mou

#Purpose: This code plots the F-statistic from PERMANOVA pairwise comparisons of INFnm to 
#either INFinject or INFfeed groups, over time (displays the magnitude of change in the fecal bacterial community structure 
#of the two antibiotic treatment groups relative to infected control)

#Files needed:
#FS9_Q2_MagnitudeOfChange.csv

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library(ggrepel)

#The F-values were obtained from FS9.WithinDayPairwiseComparisons.Q2.doubletons.txt data generated in the "Beta diversity" section
#FS9_Q2_MagnitudeOfChange.csv was created from FS9.WithinDayPairwiseComparisons.Q2.doubletons.txt data that was rearranged
#To make the FS9_Q2_MagnitudeOfChange.csv file, open FS9.WithinDayPairwiseComparisons.Q2.doubletons.txt file 
#(created from FS9_alpha_beta_diversity_Q2.R) in excel, 
#copy columns "F. Model" through "p.adjusted2" and paste in a separate spreadsheet. 
#Add "Day" and "Treatment" columns and save as "FS9_Q2_MagnitudeOfChange.csv".

fecal <- read.csv("data/FS9_Q2_MagnitudeOfChange.csv")
class(fecal)
fecal$Day <- factor(fecal$Day) #Encode "Day" as a factor
fecal$p.adjusted <- round(fecal$p.adjusted, 3)
fecal2 <- ggplot(data=fecal, aes(x=Day, y=F.Model, group=Treatment)) +
  #geom_line(aes(color=Treatment)) + 
  geom_point(aes(color=Treatment, size = 10)) +
  ylab("PERMANOVA F vs INFnm \n(difference relative to INFnm)") +
  scale_fill_manual(values = c(INFnm='#CC0066', INFinject='#E69F00', INFfeed='#999999')) +
  scale_color_manual(values = c(INFnm='#CC0066', INFinject='#E69F00', INFfeed='#999999')) +
  geom_label_repel(aes(label=p.adjusted), box.padding = 0.35, point.padding=0.5,segment.color = 'grey50', size = 5) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
  labs(color="Treatment group") +
  scale_size(guide = 'none')
fecal2

#Save 'fecal2' as a .tiff for publication, 500dpi
ggsave("FS9_Q2_Fecal_Magnitude_NoLine.tiff", plot=fecal2, width = 6, height = 5, dpi = 500, units =c("in"))
