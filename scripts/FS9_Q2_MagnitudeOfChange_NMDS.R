#####################################################################################################
#FS9 Magnitude of Change and NMDS plots combined into one figure for Q2 days 7, 11, 14 - Infected: nm vs feed or inject
#Kathy Mou

#Purpose: Combine NMDS figure and PERMANOVA F magnitude of change figure for Q2 day 7, 11, 14 into one using cowplot package

#File needed:
#FS9_Q2_MagnitudeOfChange.csv

#Load library packages
library(vegan)
library(tidyverse)
library(phyloseq)
library(philentropy)
library(cowplot)
library(ggplot2)
library(ggrepel)

#NMDS function
NMDS_ellipse <- function(metadata, OTU_table, grouping_set,
                         distance_method = 'bray',
                         rand_seed = 77777,
                         MDS_trymax = 1000,
                         autotransform = FALSE,
                         wascores = TRUE,
                         expand = FALSE){
  require(vegan)
  require(tidyr)
  
  if (grouping_set %in% colnames(metadata)){
    if (all(rownames(metadata) == rownames(OTU_table))){
      
      set.seed(rand_seed)
      generic_MDS <- metaMDS(OTU_table, k = 2,
                             trymax = MDS_trymax,
                             autotransform = autotransform,
                             distance = distance_method,
                             wascores = wascores,
                             expand = expand)
      
      stress <- generic_MDS$stress
      nmds_points <- as.data.frame(generic_MDS$points)
      metadata <- metadata[match(rownames(generic_MDS$points), rownames(metadata)),]
      metadata <- as.data.frame(metadata) # weird things were happening when a grouped tibble was used as metadata...
      metanmds <- cbind(metadata, nmds_points)
      # browser()
      nmds.mean <- aggregate(metanmds[,grep("MDS", colnames(metanmds))], list(group=metanmds[[grouping_set]]), mean)
      metanmds[[grouping_set]] <- factor(metanmds[[grouping_set]]) # this 'set' needs to be passed in from function
      
      #check to make sure at least 3 obs for each grouping_set
      
      numobs <- metanmds %>% group_by(!!grouping_set) %>% summarise(n=n())
      if (min(numobs$n) >= 3){
        ord <- ordiellipse(generic_MDS, metanmds[[grouping_set]], label = TRUE, conf = .95, kind = 'se', draw = 'none')
        
        df_ell <- data.frame()
        for (d in levels(metanmds[[grouping_set]])){
          df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds[[grouping_set]] == d,],
                                                           vegan:::veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
        }
        
        # this loop assigns the group centroid X coordinates to each sample, there is probably a better way...
        
        metanmds$centroidX <- NA
        metanmds$centroidY <- NA
        
        for (level in levels(metanmds[[grouping_set]])){
          metanmds[metanmds[[grouping_set]] == level,]$centroidX <- nmds.mean$MDS1[nmds.mean$group == level]
          metanmds[metanmds[[grouping_set]] == level,]$centroidY <- nmds.mean$MDS2[nmds.mean$group == level]
          
        }
        print(paste('Ordination stress:', stress, sep = ' '))
        return(list(metanmds, df_ell, generic_MDS))
        
      } else {
        warning('One of your groups in "grouping_set" has less than 3 observations, cannot generate elipses')
        df_ell <- data.frame()
        return(list(metanmds, df_ell, generic_MDS))}
      
    } else {
      stop('The rownames for your OTU table and metadata do not match.')
    }
    
  }else {
    stop('The grouping set column you have provided in not in your metadata.')
  }
}

################################################
#Run FS9_phyloseq_Q2.R to generate phyloseq.FS9 object to run the R script below.

#Setting up 'phyloseq' into dataframes for NMDS calculation
meta <- data.frame(phyloseq.FS9@sam_data) #Make 'phyloseq.FS9' sam_data into dataframe
otu <- data.frame((phyloseq.FS9@otu_table)) #Make 'phyloseq.FS9' otu_table into dataframe
class(meta) #data.frame
rownames(meta) == row.names(otu) #Make sure rownames between 'meta' and 'otu' match exactly. It is true
meta$numOTUS <- rowSums(otu > 1) #For rows with sums greater than 1 in 'otu', move rows and their respective sum values into "numOTUs" column in 'meta'
head(meta)

#NMDS calculation (aka beta diversity)
otu[1:10,1:10]
dim(otu) #72 1042 (doubletons removed)
NMDS <- NMDS_ellipse(meta, otu, grouping_set = 'All')
#Output:
#[1] "Ordination stress: 0.181970397207524"

#Separate meta data and ellipse data to two lists to make NMDS plot
head(NMDS)
metanmds <- NMDS[[1]] #'metanmds' has meta data + MDS calculations. Select this 1st list of 'NMDS' using double brackets
df_ell <- NMDS[[2]] #'df_ell' is accessing 2nd list from 'NMDS' that has ellipse calculations

#Need two separate lists for making NMDS plot
df_ell$group
head(df_ell)

#Create "Day" and "Treatment" columns within 'df_ell' for faceting purposes
df_ell <- df_ell %>% separate(group, into=c("Day","Treatment"), sep="_", remove=FALSE)
View(df_ell)

#Restructure level order for 'metanmds' and 'df_ell'
unique(metanmds$Day) #"D11" "D7"  "D14"
unique(df_ell$Day) #"D11" "D7"  "D14"
metanmds$Day = factor(metanmds$Day, levels = c("D7", "D11", "D14"))
df_ell$Day = factor(df_ell$Day, levels = c("D7", "D11", "D14"))
levels(df_ell$Day) #"D7"  "D11" "D14"
levels(metanmds$Day) #"D7"  "D11" "D14"

#All points, all days, all treatments, samples that aren't relevant are grayed out on plot
#All points, all days, all treatments, samples that aren't relevant are grayed out on plot
#Plotting with gridlines and axes, gray points for All days
metanmds.2 <- metanmds
metanmds.2$Treatment2 = metanmds.2$Treatment
metanmds.2$Treatment2 <- as.character(metanmds.2$Treatment2)
metanmds.2$Treatment2

#All days and treatments faceted by day (gridlines)
nmdsplot_treatment2<- ggplot(metanmds, aes(x=MDS1, y=MDS2)) +  annotate(x=metanmds.2$MDS1, y=metanmds.2$MDS2, color='grey57', geom = 'point')+
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment), size=1.25) + 
  geom_point(aes(color = Treatment), size = 2) + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=Treatment), alpha=.5) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + facet_wrap(~Day) +
  theme_bw() +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  labs(caption = 'Ordination stress = 0.182', color="Treatment group")
nmdsplot_treatment2

#Make magnitude of change plot

#Purpose: This code plots the F-statistic from PERMANOVA pairwise comparisons of INFnm to 
#either INFinject or INFfeed groups, over time (displays the magnitude of change in the fecal bacterial community structure 
#of the two antibiotic treatment groups relative to infected control)

#Files needed:
#FS9_Q2_MagnitudeOfChange.csv

#The F-values were obtained from FS9.WithinDayPairwiseComparisons.Q2.doubletons.txt data generated in the "Beta diversity" section
#FS9_Q2_MagnitudeOfChange.csv was created from FS9.WithinDayPairwiseComparisons.Q2.doubletons.txt data that was rearranged
#To make the FS9_Q2_MagnitudeOfChange.csv file, open FS9.WithinDayPairwiseComparisons.Q2.doubletons.txt file 
#(created from FS9_alpha_beta_diversity_Q2.R) in excel, 
#copy columns "F. Model" through "p.adjusted2" and paste in a separate spreadsheet. 
#Add "Day" and "Treatment" columns and save as "FS9_Q2_MagnitudeOfChange.csv".

fecal <- read.csv("./data/FS9_Q2_MagnitudeOfChange.csv")
class(fecal)
fecal$Day <- factor(fecal$Day) #Encode "Day" as a factor
fecal$p.adjusted <- round(fecal$p.adjusted, 3)
fecal2 <- ggplot(data=fecal, aes(x=Day, y=F.Model, group=Treatment)) +
  #geom_line(aes(color=Treatment)) + 
  geom_point(aes(color=Treatment, size = 10)) +
  ylab("PERMANOVA F vs INFnm") +
  scale_fill_manual(values = c(INFinject='#00BA38',
                               INFnm='#F8766D',
                               INFfeed='#619CFF')) +
  scale_color_manual(values = c(INFinject='#00BA38',
                                INFnm='#F8766D',
                                INFfeed='#619CFF')) +
  geom_label_repel(aes(label=p.adjusted), box.padding = 0.35, point.padding=0.5,segment.color = 'grey50', size = 5) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size=16), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), legend.text=element_text(size=16), legend.title=element_text(size=16)) +
  labs(color="Treatment group") +
  scale_size(guide = 'none') +
  theme_bw()
fecal2

#Cowplot of Q1_NONINFnm_INFnm Order and Phylum DNEG3, D7
Fig6ab <- plot_grid(fecal2, nmdsplot_treatment2, labels = c('A', 'B'), label_size = 12, ncol=1)
#, rel_widths = c(1, 4)
Fig6ab

ggsave("Q2_INFnm_inject_feed_MagnitudeOfChange_NMDS.tiff", plot=Fig6ab, width = 10, height = 9, dpi = 500, units =c("in"))
