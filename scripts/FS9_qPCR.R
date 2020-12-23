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

tet32stats <- tet32 %>% group_by(Day,Treatment) %>% summarise(mean=mean(value), 
                                                                   n=n(), 
                                                                   sd=sd(value), 
                                                                   se=sd/n) %>% write_csv('./results/tet32_mean.csv')

tet32 %>% 
  ggplot(aes(x=Day, y=log10tet32, color=Treatment)) +
  scale_color_manual(values = c(INFinject='#E69F00', INFnm='#CC0066', INFfeed='#999999')) +
  geom_boxplot() + 
  ylab('log10-fold change of tet32 gene abundance') +
  #axis label font size! and axis tick font size!
  theme_bw()
