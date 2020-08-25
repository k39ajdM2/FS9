#######################################################################
#FS9 16S alpha and beta diversity - Noninfected vs Infected, additional analysis
#Kathy Mou

#NOTES: 
#This code analyzes alpha and beta diversity statistics for fecal samples from all days or only days 4 and 7; and associated plots
#This script uses files created in "FS9_phyloseq.R"

#Clear workspace and load necessary packages
rm(list=ls())

#Load library packages
library(vegan)
library(tidyverse)
library(phyloseq)
library(scales)
library(RColorBrewer)
library(philentropy)
library(cowplot)

#Load NMDS_ellipse function from https://github.com/Jtrachsel/funfuns
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


#Load pairwise.adonis function from https://github.com/Jtrachsel/funfuns, I checked with Jules' code on 1May2020 and it is up-to-date to what he has in funfuns
#this function taken from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m = 'none', permutations=999){
  require(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = permutations);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

###########################################################################################################
#Load 'phyloseq.FS9.RData' into environment
load('./data/phyloseq.FS9.doubleton.RData')

#Setting up 'phyloseq' into dataframes for NMDS calculation
meta <- data.frame(phyloseq.FS9@sam_data) #Make 'phyloseq.FS9' sam_data into dataframe
otu <- data.frame((phyloseq.FS9@otu_table)) #Make 'phyloseq.FS9' otu_table into dataframe
class(meta) #data.frame
rownames(meta) == row.names(otu) #Make sure rownames between 'meta' and 'otu' match exactly. It is true
meta$numOTUS <- rowSums(otu > 1) #For rows with sums greater than 1 in 'otu', move rows and their respective sum values into "numOTUs" column in 'meta'
head(meta)

#NMDS calculation (aka beta diversity)
otu[1:10,1:10]
dim(otu) #79 1156 (doubletons removed)
NMDS <- NMDS_ellipse(meta, otu, grouping_set = 'All')
#Output:
#[1] "Ordination stress: 0.177167945792148"

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
unique(metanmds$Day) #"D4"    "DNEG3" "D0"    "D7" 
unique(df_ell$Day) #"D0"    "D4"    "D7"    "DNEG3" 
metanmds$Day = factor(metanmds$Day, levels = c("DNEG3", "D0", "D4", "D7"))
df_ell$Day = factor(df_ell$Day, levels = c("DNEG3", "D0", "D4", "D7"))
levels(df_ell$Day) #"DNEG3" "D0"    "D4"    "D7"
levels(metanmds$Day) #"DNEG3" "D0"    "D4"    "D7"
dim(metanmds) #79 10
dim(df_ell) #808 5

#Creating NMDS day+treatment plot from NMDS calculations
nmdsplot <- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) + 
  facet_wrap(~Day, scales = 'free') +
  #scale_color_brewer(palette="Dark2") +
  theme_gray(base_size = 10) +
  theme(strip.text.x = element_text(size=15), axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
  labs(color="Treatment group")+
  scale_color_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  labs(caption = 'Ordination stress = 0.177')
#nmdsplot2 <- nmdsplot + scale_colour_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position = "right")
nmdsplot
#Save 'nmdsplot' as a .tiff for publication, 500dpi
#ggsave("NMDS_DayAndTreatment.tiff", plot=nmdsplot, width = 11, height = 5, dpi = 500, units =c("in"))

#Creating NMDS day plot from NMDS calculations
nmdsplot_day <- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Day)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = .5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Day, group=group)) + 
  labs(caption = 'Ordination stress = 0.177') 
nmdsplot_day
#Save 'nmdsplot_day' as a .tiff for publication, 500dpi
#ggsave("NMDS_Day.tiff", plot=nmdsplot_day, width = 10, height = 6, dpi = 500, units =c("in"))

#Creating NMDS treatment plot from NMDS calculations
nmdsplot_treatment<- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = .5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) + 
  scale_color_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  labs(caption = 'Ordination stress = 0.177')  
nmdsplot_treatment
#Save 'nmdsplot_treatment' as a .tiff for publication, 500dpi
#ggsave("NMDS_Treatment.tiff", plot=nmdsplot_treatment, width = 10, height = 6, dpi = 500, units =c("in"))

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
  scale_color_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  labs(caption = 'Ordination stress = 0.177', color="Treatment group")
nmdsplot_treatment2
#Save 'nmdsplot_treatment2' as a .tiff for publication, 500dpi
#ggsave("NMDS_DayAndTreatment_AllSamples.tiff", plot=nmdsplot_treatment2, width = 10, height = 6, dpi = 500, units =c("in"))

#Using pairwise.adonis function (beta diversity)
adon <- pairwise.adonis(otu, meta$All) #Run pairwise.adonis on 'otu' OTU table and "All" column of 'meta' dataframe
#adon contains all the pairwise comparisons
adon$pairs #List all comparisons in the "pairs" column of 'nw.adon'
goodcomps <- c(grep('DNEG3_[A-Za-z]+ vs DNEG3_[A-Za-z]+', adon$pairs),
               grep('D0_[A-Za-z]+ vs D0_[A-Za-z]+', adon$pairs),
               grep('D4_[A-Za-z]+ vs D4_[A-Za-z]+', adon$pairs),
               grep('D7_[A-Za-z]+ vs D7_[A-Za-z]+', adon$pairs))
# "[A-Za-z]" matches all capital and lowercase letters
# "+" matches a whole word and not just one letter (if you didn't have "+", then it would match by one letter)
# "c" creates the vector, lumps all pairs of specific groups of interest together
# You want to make a vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
adon.good <- adon[goodcomps,] #Rename 'goodcomps' vector to 'adon.good'
adon.good
adon.good$p.adjusted <- p.adjust(adon.good$p.value, method = 'fdr') #"p.adjust" function returns a set of p-values adjusted with "fdr" method
adon.good$p.adjusted2 <- round(adon.good$p.adjusted, 3) #Round p-values to 3 decimal points and list in new "p.adjusted2" column
adon.good$p.adjusted2[adon.good$p.adjusted2 > 0.05] <- NA #For all p-values greater than 0.05, replace with "NA"
write.csv(adon.good, file='FS9.WithinDayPairwiseComparisons.Q1.doubletons.alldays.txt', row.names=TRUE)

#Alpha diversity
#Calculating alpha diversity metrics: Shannon, Inverse Simpson
meta$shannon <- diversity(otu) #"diversity" is a vegan function. The default index is set at "shannon". I added a shannon index column in 'meta'
meta$invsimpson <- diversity(otu,index = 'invsimpson') #We used 'invsimpson' since it is easier to interpret than Simpson values and won't need to "inverse" the Simpson values to understand (With Simpson values, the lower the number, the higher the diversity)
levels(sample_data(meta)$Day)
meta$Day = factor(meta$Day, levels = c("DNEG3", "D0", "D4", "D7"))  # Set the level order of values in "Day" column
levels(sample_data(meta)$Day) #"DNEG3" "D0"  "D4"    "D7"  

#Calculate the average shannon, invsimpson, numOTUs for each "All" subtype within meta
shannon.invsimpson.numOTUs <- aggregate(meta[, 6:8], list(meta$All), mean)
print(shannon.invsimpson.numOTUs)
#Output (singletons removed):
#        Group.1    numOTUS   shannon    invsimpson
#1       D0_INFnm  83.00000 3.153554  15.628239
2    D0_NONINFnm  71.33333 3.089266  13.274195
3       D4_INFnm  79.66667 3.202696  11.903607
4    D4_NONINFnm  77.60000 3.117897   9.115051
5       D7_INFnm  94.80000 3.661901  19.286844
6    D7_NONINFnm  92.37500 3.770617  18.321736
7    DNEG3_INFnm 113.11765 4.062337  31.949254
8 DNEG3_NONINFnm  97.64706 3.864239  25.325367
write.csv(shannon.invsimpson.numOTUs, file="FS9.shannon.invsimpson.num.OTUs.doubletons.txt", row.names=TRUE)

#Shannon
pairwise.wilcox.shannon.test <- pairwise.wilcox.test(meta$shannon, meta$All, p.adjust.method = 'none') #Calculate pairwise comparisons by "All" column of the shannon indices in "Shannon" column
print(pairwise.wilcox.shannon.test) #Look at the results of 'pairwise.wilcox.shannon.test'

#Inverse Simpson
pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(meta$invsimpson, meta$All, p.adjust.method = 'none')
print(pairwise.wilcox.invsimpson.test)

#Generate a box and whisker plot of shannon (both shannon and inverse simpson diversity indices showed same trends: 
#no significant differences between treatment groups within a day)
shan <- ggplot(data = meta, aes(x=All, y=shannon, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Shannon diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text.x = element_text(size=14),
        axis.title.y = element_text(size=15)) +
  scale_fill_manual(values = c(INF='#CC0066', NONINF='#56B4E9')) +
  theme(legend.position = "none")
shan
# "free" within "facet_wrap" allows each plot to customize the scale to the specific data set (no forced scaling applied to all plots)
# "position = position_dodge2(preserve = 'total')" fixes the ggplot box width, making them wider, prevents narrow boxes from forming in the plot

#Save 'shan' as a .tiff for publication, 500dpi
#ggsave("FS9_Shannon.tiff", plot=shan, width = 7, height = 7, dpi = 500, units =c("in"))

#Generate a box and whisker plot of inverse simpson 
invsimp <- ggplot(data = meta, aes(x=All, y=invsimpson, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Inverse Simpson diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text.x = element_text(size=14),
        axis.title.y = element_text(size=15)) +
  scale_fill_manual(values = c(INF='#CC0066', NONINF='#56B4E9')) +
  theme(legend.position = "none")
invsimp
#Save 'invsimp' as a .tiff for publication, 500dpi
#ggsave("FS9_InverseSimpson.tiff", plot=invsimp, width = 7, height = 7, dpi = 500, units =c("in"))

alphadiv <- plot_grid(shan, invsimp, labels = c('A', 'B'), label_size = 12) #How to fix plot dimensions so that y-axis labels doesn't overlap the plots?
ggsave(alphadiv,
       filename = './results/alpha_diversity_NONINFvsINF.jpeg',
       width = 300,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')





###################################################################### Days 4, 7 only ##########
#Setting up 'phyloseq' into dataframes for NMDS calculation
meta <- data.frame(phyloseq.FS9@sam_data) #Make 'phyloseq.FS9' sam_data into dataframe
otu <- data.frame((phyloseq.FS9@otu_table)) #Make 'phyloseq.FS9' otu_table into dataframe
class(meta) #data.frame
rownames(meta) == row.names(otu) #Make sure rownames between 'meta' and 'otu' match exactly. It is true
meta$numOTUS <- rowSums(otu > 1) #For rows with sums greater than 1 in 'otu', move rows and their respective sum values into "numOTUs" column in 'meta'
head(meta)

#NMDS calculation (aka beta diversity)
otu[1:10,1:10]
dim(otu) #29 768 (doubletons removed)
NMDS <- NMDS_ellipse(meta, otu, grouping_set = 'All')
#Output:
#[1] "Ordination stress: 0.171192923672683"

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
unique(metanmds$Day) #"D4" "D7" 
unique(df_ell$Day) #""D4"    "D7"
metanmds$Day = factor(metanmds$Day, levels = c("D4", "D7"))
df_ell$Day = factor(df_ell$Day, levels = c("D4", "D7"))
levels(df_ell$Day) #"D4" "D7"
levels(metanmds$Day) #"D4" "D7"
dim(metanmds) #29 10
dim(df_ell) #404 5

#Creating NMDS day+treatment plot from NMDS calculations
nmdsplot <- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) + 
  facet_wrap(~Day, scales = 'free') +
  #scale_color_brewer(palette="Dark2") +
  theme_gray(base_size = 10) +
  theme(strip.text.x = element_text(size=15), axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
  labs(color="Treatment group")+
  scale_color_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  labs(caption = 'Ordination stress = 0.171')
#nmdsplot2 <- nmdsplot + scale_colour_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position = "right")
nmdsplot
#Save 'nmdsplot' as a .tiff for publication, 500dpi
#ggsave("NMDS_DayAndTreatment.tiff", plot=nmdsplot, width = 11, height = 5, dpi = 500, units =c("in"))

#Creating NMDS day plot from NMDS calculations
nmdsplot_day <- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Day)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = .5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Day, group=group)) + 
  labs(caption = 'Ordination stress = 0.171') 
nmdsplot_day
#Save 'nmdsplot_day' as a .tiff for publication, 500dpi
#ggsave("NMDS_Day.tiff", plot=nmdsplot_day, width = 10, height = 6, dpi = 500, units =c("in"))

#Creating NMDS treatment plot from NMDS calculations
nmdsplot_treatment<- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = .5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) + 
  scale_color_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  labs(caption = 'Ordination stress = 0.171')  
nmdsplot_treatment
#Save 'nmdsplot_treatment' as a .tiff for publication, 500dpi
#ggsave("NMDS_Treatment.tiff", plot=nmdsplot_treatment, width = 10, height = 6, dpi = 500, units =c("in"))

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
  scale_color_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  labs(caption = 'Ordination stress = 0.171', color="Treatment group")
nmdsplot_treatment2
#Save 'nmdsplot_treatment2' as a .tiff for publication, 500dpi
#ggsave("NMDS_DayAndTreatment_AllSamples.tiff", plot=nmdsplot_treatment2, width = 10, height = 6, dpi = 500, units =c("in"))

#Using pairwise.adonis function (beta diversity)
adon <- pairwise.adonis(otu, meta$All) #Run pairwise.adonis on 'otu' OTU table and "All" column of 'meta' dataframe
#adon contains all the pairwise comparisons
adon$pairs #List all comparisons in the "pairs" column of 'nw.adon'
goodcomps <- c(grep('D4_[A-Za-z]+ vs D4_[A-Za-z]+', adon$pairs),
               grep('D7_[A-Za-z]+ vs D7_[A-Za-z]+', adon$pairs))
# "[A-Za-z]" matches all capital and lowercase letters
# "+" matches a whole word and not just one letter (if you didn't have "+", then it would match by one letter)
# "c" creates the vector, lumps all pairs of specific groups of interest together
# You want to make a vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
adon.good <- adon[goodcomps,] #Rename 'goodcomps' vector to 'adon.good'
adon.good
adon.good$p.adjusted <- p.adjust(adon.good$p.value, method = 'fdr') #"p.adjust" function returns a set of p-values adjusted with "fdr" method
adon.good$p.adjusted2 <- round(adon.good$p.adjusted, 3) #Round p-values to 3 decimal points and list in new "p.adjusted2" column
adon.good$p.adjusted2[adon.good$p.adjusted2 > 0.05] <- NA #For all p-values greater than 0.05, replace with "NA"
write.csv(adon.good, file='FS9.WithinDayPairwiseComparisons.Q1.doubletons.D4D7.txt', row.names=TRUE)

#Alpha diversity
#Calculating alpha diversity metrics: Shannon, Inverse Simpson
meta$shannon <- diversity(otu) #"diversity" is a vegan function. The default index is set at "shannon". I added a shannon index column in 'meta'
meta$invsimpson <- diversity(otu,index = 'invsimpson') #We used 'invsimpson' since it is easier to interpret than Simpson values and won't need to "inverse" the Simpson values to understand (With Simpson values, the lower the number, the higher the diversity)
levels(sample_data(meta)$Day)
meta$Day = factor(meta$Day, levels = c("D4", "D7"))  # Set the level order of values in "Day" column
levels(sample_data(meta)$Day) #"D4" "D7"

#Calculate the average shannon, invsimpson, numOTUs for each "All" subtype within meta
shannon.invsimpson.numOTUs <- aggregate(meta[, 6:8], list(meta$All), mean)
print(shannon.invsimpson.numOTUs)
#Output (singletons removed):
#        Group.1    numOTUS   shannon    invsimpson
#1    D4_INFnm    79.66667    3.202696  11.903607
#2 D4_NONINFnm    77.60000    3.117897   9.115051
#3    D7_INFnm    94.80000    3.661901  19.286844
#4 D7_NONINFnm    92.37500    3.770617  18.321736
write.csv(shannon.invsimpson.numOTUs, file="FS9.shannon.invsimpson.num.OTUs.doubletons.txt", row.names=TRUE)

#Shannon
pairwise.wilcox.shannon.test <- pairwise.wilcox.test(meta$shannon, meta$All, p.adjust.method = 'none') #Calculate pairwise comparisons by "All" column of the shannon indices in "Shannon" column
print(pairwise.wilcox.shannon.test) #Look at the results of 'pairwise.wilcox.shannon.test'

#Inverse Simpson
pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(meta$invsimpson, meta$All, p.adjust.method = 'none')
print(pairwise.wilcox.invsimpson.test)

#Generate a box and whisker plot of shannon (both shannon and inverse simpson diversity indices showed same trends: 
#no significant differences between treatment groups within a day)
shan <- ggplot(data = meta, aes(x=All, y=shannon, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Shannon diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text.x = element_text(size=14),
        axis.title.y = element_text(size=15)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  theme(legend.position = "none")
shan
# "free" within "facet_wrap" allows each plot to customize the scale to the specific data set (no forced scaling applied to all plots)
# "position = position_dodge2(preserve = 'total')" fixes the ggplot box width, making them wider, prevents narrow boxes from forming in the plot

#Save 'shan' as a .tiff for publication, 500dpi
#ggsave("FS9_Shannon.tiff", plot=shan, width = 7, height = 7, dpi = 500, units =c("in"))

#Generate a box and whisker plot of inverse simpson 
invsimp <- ggplot(data = meta, aes(x=All, y=invsimpson, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Inverse Simpson diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text.x = element_text(size=14),
        axis.title.y = element_text(size=15)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  theme(legend.position = "none")
invsimp
#Save 'invsimp' as a .tiff for publication, 500dpi
#ggsave("FS9_InverseSimpson.tiff", plot=invsimp, width = 7, height = 7, dpi = 500, units =c("in"))

alphadiv <- plot_grid(shan, invsimp, labels = c('A', 'B'), label_size = 12) #How to fix plot dimensions so that y-axis labels doesn't overlap the plots?
ggsave(alphadiv,
       filename = './results/alpha_diversity_NONINFvsINF_D4D7.jpeg',
       width = 300,
       height = 120,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')
