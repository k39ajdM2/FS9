#######################################################################
#FS9 16S alpha and beta diversity
#Kathy Mou

#NOTES: 
#This code analyzes alpha and beta diversity statistics for fecal samples, and associated plots
#This script uses files created in "FS9_phyloseq.R"

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory (either on desktop or network drive)
setwd("~/Desktop/FS9/FS9_RWorkspace")


#Load library packages
install.packages("vegan")
install.packages("tidyverse")
install.packages("phyloseq")
install.packages("scales")
install.packages("RColorBrewer")
install.packages("philentropy")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(vegan)
library(tidyverse)
library(phyloseq)
library(scales)
library(RColorBrewer)
library(philentropy)

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

#Load image file
load("FS9_alpha_beta_diversity.RData")

#Save image file
save.image(file="FS9_alpha_beta_diversity.RData")

###########################################################################################################
#Load 'FS9.phyloseq.RData' into environment by clicking on this file name in the correct directory of the Files/Plots/Packages/Help panel

#Setting up 'phyloseq' into dataframes for NMDS calculation
meta <- data.frame(phyloseq@sam_data) #Make 'phyloseq' sam_data into dataframe
otu <- data.frame(t(phyloseq@otu_table)) #Make 'phyloseq' otu_table into dataframe
class(meta) #data.frame
rownames(meta) == row.names(otu) #Make sure rownames between 'meta' and 'otu' match exactly. It is true
meta$numOTUS <- rowSums(otu > 1) #For rows with sums greater than 1 in 'otu', move rows and their respective sum values into "numOTUs" column in 'meta'
head(meta)

#NMDS calculation
otu[1:10,1:10]
NMDS <- NMDS_ellipse(meta, otu, grouping_set = 'All')
#Output:
#Result: [1] "Ordination stress: 0.198470463533747"

#Separate meta data and ellipse data to two lists to make NMDS plot
head(NMDS)
metanmds <- NMDS[[1]] #'metanmds' has meta data + MDS calculations. Select this 1st list of 'NMDS' using double brackets
df_ell <- NMDS[[2]] #'df_ell' is accessing 2nd list from 'NMDS' that has ellipse calculations

#Need two separate lists for making NMDS plot
df_ell$group
head(df_ell)

#Create "Day" and "Treatment" columns within 'df_ell' for faceting purposes
df_ell <- df_ell %>% separate(group, into=c("Day","Treatment"), sep=" ", remove=FALSE) #Daniel's good stuff
View(df_ell)

#Restructure level order for 'metanmds' and 'df_ell'
unique(metanmds$Day) #-3 0 7
unique(df_ell$Day) #"-3" "0"  "7" 
metanmds$Day = factor(metanmds$Day, levels = c("-3", "0", "7"))
df_ell$Day = factor(df_ell$Day, levels = c("-3", "0", "7"))
levels(df_ell$Day) #"-3" "0"  "7" 
levels(metanmds$Day) #"-3" "0"  "7" 

#Skip renaming treatment groups until have more info
#Renaming treatment groups xxx to xxx, respectively, in 'metanmds' and 'df_ell' dataframes
metanmds$Treatment2 = metanmds$Treatment
metanmds$Treatment2 <- as.character(metanmds$Treatment2)
metanmds$Treatment2[metanmds$Treatment2 == 'INF_OralOCT'] <- "INF_OralOTC"
metanmds$Treatment2[metanmds$Treatment2 == 'INF_InjOCT'] <- "INF_InjOTC"
df_ell$Treatment2 = df_ell$Treatment
df_ell$Treatment2 <- as.character(df_ell$Treatment2)
df_ell$Treatment2[df_ell$Treatment2 == 'INF_OralOCT'] <- "INF_OralOTC"
df_ell$Treatment2[df_ell$Treatment2 == 'INF_InjOCT'] <- "INF_InjOTC"
metanmds$Treatment2
df_ell$Treatment2

#Creating NMDS day+treatment plot from NMDS calculations
nmdsplot <- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Treatment2)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment2, group=group)) + 
  facet_wrap(~Day, scales = 'free') +
  #scale_color_brewer(palette="Dark2") +
  theme_gray(base_size = 10) +
  theme(strip.text.x = element_text(size=15), axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
  labs(color="Treatment group")+
  labs(caption = 'Ordination stress = 0.198')
#nmdsplot2 <- nmdsplot + scale_colour_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position = "right")
nmdsplot
#Save 'nmdsplot' as a .tiff for publication, 500dpi
ggsave("NMDS_DayAndTreatment.tiff", plot=nmdsplot, width = 11, height = 5, dpi = 500, units =c("in"))


#Creating NMDS day plot from NMDS calculations
nmdsplot_day <- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Day)) + geom_point() + 
    geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = .5) + 
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Day, group=group)) + 
    labs(caption = 'Ordination stress = 0.198') 
nmdsplot_day
#Save 'nmdsplot_day' as a .tiff for publication, 500dpi
ggsave("NMDS_Day.tiff", plot=nmdsplot_day, width = 10, height = 6, dpi = 500, units =c("in"))


#Creating NMDS treatment plot from NMDS calculations
nmdsplot_treatment<- ggplot(data=metanmds, aes(x=MDS1, y=MDS2, color=Treatment2)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = .5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment2, group=group)) + 
  labs(caption = 'Ordination stress = 0.198')  
nmdsplot_treatment
#Save 'nmdsplot_treatment' as a .tiff for publication, 500dpi
ggsave("NMDS_Treatment.tiff", plot=nmdsplot_treatment, width = 10, height = 6, dpi = 500, units =c("in"))


#All points, all days, all treatments, gray points that aren't relevant

#Plotting with gridlines and axes, gray points for All days
metanmds.2 <- metanmds
metanmds.2$Treatment2 = metanmds.2$Treatment
metanmds.2$Treatment2 <- as.character(metanmds.2$Treatment2)
metanmds.2$Treatment2[metanmds.2$Treatment2 == 'INF_OralOCT'] <- "INF_OralOTC"
metanmds.2$Treatment2[metanmds.2$Treatment2 == 'INF_InjOCT'] <- "INF_InjOTC"
metanmds.2$Treatment2


#All days and treatments faceted by day (gridlines)
nmdsplot_treatment2<- ggplot(metanmds, aes(x=MDS1, y=MDS2)) +  annotate(x=metanmds.2$MDS1, y=metanmds.2$MDS2, color='grey57', geom = 'point')+
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment2), size=1.25) + 
  geom_point(aes(color = Treatment2), size = 2) + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=Treatment2), alpha=.5) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + facet_wrap(~Day, nrow = 1) +
  theme_bw() +
  labs(caption = 'Ordination stress = 0.198', color="Treatment group")
nmdsplot_treatment2
#Save 'nmdsplot_treatment2' as a .tiff for publication, 500dpi
ggsave("NMDS_DayAndTreatment_AllSamples.tiff", plot=nmdsplot_treatment2, width = 10, height = 6, dpi = 500, units =c("in"))


#Using pairwise.adonis function
adon <- pairwise.adonis(otu, meta$All) #Run pairwise.adonis on 'otu' OTU table and "All" column of 'meta' dataframe
#adon contains all the pairwise comparisons
adon$pairs #List all comparisons in the "pairs" column of 'nw.adon'
goodcomps <- c(grep('-3 [A-Za-z]+ vs -3 [A-Za-z]+', adon$pairs),
                  grep('0 [A-Za-z]+ vs 0 [A-Za-z]+', adon$pairs),
                  grep('7 [A-Za-z]+ vs 7 [A-Za-z]+', adon$pairs))
# "[A-Za-z]" matches all capital and lowercase letters
# "+" matches a whole word and not just one letter (if you didn't have "+", then it would match by one letter)
# "c" creates the vector, lumps all pairs of specific groups of interest together
# You want to make a vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
adon.good <- adon[goodcomps,] #Rename 'goodcomps' vector to 'adon.good'
adon.good
adon.good$p.adjusted <- p.adjust(adon.good$p.value, method = 'fdr') #"p.adjust" function returns a set of p-values adjusted with "fdr" method
adon.good$p.adjusted2 <- round(adon.good$p.adjusted, 3) #Round p-values to 3 decimal points and list in new "p.adjusted2" column
adon.good$p.adjusted2[adon.good$p.adjusted2 > 0.05] <- NA #For all p-values greater than 0.05, replace with "NA"
write.csv(adon.good, file='FS8b.adon.good.txt', row.names=TRUE)
#4 pairs of significant differences in microbial composition

adon$p.adjusted <- p.adjust(adon$p.value, method = 'fdr') #"p.adjust" function returns a set of p-values adjusted with "fdr" method
adon$p.adjusted2 <- round(adon$p.adjusted, 3) #Round p-values to 3 decimal points and list in new "p.adjusted2" column
adon$p.adjusted2[adon$p.adjusted2 > 0.05] <- NA #For all p-values greater than 0.05, replace with "NA"
write.csv(adon, file='FS9.adon.txt', row.names=TRUE)


#Alpha diversity
#Calculating alpha diversity metrics: Shannon, Inverse Simpson
meta$shannon <- diversity(otu) #"diversity" is a vegan function. The default index is set at "shannon". I added a shannon index column in 'meta'
meta$invsimpson <- diversity(otu,index = 'invsimpson') #We used 'invsimpson' since it is easier to interpret than Simpson values and won't need to "inverse" the Simpson values to understand (With Simpson values, the lower the number, the higher the diversity)
levels(sample_data(meta)$Day) # Set the level order of values in "Day" column

#Calculate the average shannon, invsimpson, numOTUs for each "All" subtype within meta
shannon.invsimpson.numOTUs <- aggregate(meta[, 6:8], list(meta$All), mean)
print(shannon.invsimpson.numOTUs)
#Output:
#     Group.1           numOTUS   shannon     invsimpson
#1    -3 INF_InjOTC     95.37500  3.858225    23.77920
#2    -3 INF_NoTRMT     110.62500 4.097624    33.84238
#3    -3 INF_OralOTC    92.75000  3.805127    23.79197
#4    -3 NONINF_NoTRMT  101.38889 3.930916    26.24336
#5    0 INF_InjOTC      84.23077  3.412381    16.90210
#6    0 INF_NoTRMT      84.85714  3.157894    15.53953
#7    0 INF_OralOTC     89.40000  3.540890    15.39673
#8    0 NONINF_NoTRMT   72.33333  3.060047    13.08113
#9    7 INF_InjOTC      99.30000  3.874105    22.20230
#10   7 INF_NoTRMT      95.20000  3.699163    18.93731
#11   7 INF_OralOTC     89.40000  3.628494    17.16434
#12   7 NONINF_NoTRMT   97.62500  3.785084    18.41787

write.csv(shannon.invsimpson.numOTUs, file="FS9.shannon.invsimpson.num.OTUs.txt", row.names=TRUE)
pairwise.wilcox.shannon.test <- pairwise.wilcox.test(meta$shannon, meta$All, p.adjust.method = 'none') #Calculate pairwise comparisons by "All" column of the shannon indices in "Shannon" column
print(pairwise.wilcox.shannon.test) #Look at the results of 'pairwise.wilcox.shannon.test'

#Pairwise comparisons using Wilcoxon rank sum test 

#data:  meta$shannon and meta$All 

#-3 INF_InjOTC -3 INF_NoTRMT -3 INF_OralOTC -3 NONINF_NoTRMT 0 INF_InjOTC 0 INF_NoTRMT
#-3 INF_NoTRMT    0.0938        -             -              -                -            -           
#-3 INF_OralOTC   0.8381        0.0575        -              -                -            -           
#-3 NONINF_NoTRMT 0.6212        0.1643        0.3613         -                -            -           
# 0 INF_InjOTC     0.0916        0.0035        0.1813         0.0221           -            -           
# 0 INF_NoTRMT     0.0394        0.0329        0.0414         0.0292           0.3114       -           
# 0 INF_OralOTC    0.1530        0.0318        0.2718         0.0796           0.9241       0.2677      
# 0 NONINF_NoTRMT  0.0116        0.0011        0.0300         0.0043           0.3933       0.8371      
# 7 INF_InjOTC     0.7760        0.0467        0.8121         0.3318           0.2080       0.0250      
# 7 INF_NoTRMT     0.4518        0.0684        0.5884         0.4357           0.3434       0.1613      
# 7 INF_OralOTC    0.2199        0.0087        0.3504         0.0452           0.8793       0.1088      
# 7 NONINF_NoTRMT  0.5283        0.0448        0.7461         0.1961           0.4137       0.0541      
# 0 INF_OralOTC 0 NONINF_NoTRMT 7 INF_InjOTC 7 INF_NoTRMT 7 INF_OralOTC
# -3 INF_NoTRMT    -             -               -            -            -            
# -3 INF_OralOTC   -             -               -            -            -            
# -3 NONINF_NoTRMT -             -               -            -            -            
# 0 INF_InjOTC     -             -               -            -            -            
# 0 INF_NoTRMT     -             -               -            -            -            
# 0 INF_OralOTC    -             -               -            -            -            
# 0 NONINF_NoTRMT  0.6064        -               -            -            -            
# 7 INF_InjOTC     0.1292        0.0653          -            -            -            
# 7 INF_NoTRMT     0.5135        0.0789          0.9118       -            -            
# 7 INF_OralOTC    0.7679        0.3154          0.1655       0.6842       -            
# 7 NONINF_NoTRMT  0.3543        0.0745          0.9654       0.9654       0.3599       

#P value adjustment method: none 


pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(meta$invsimpson, meta$All, p.adjust.method = 'none')
print(pairwise.wilcox.invsimpson.test)

#Pairwise comparisons using Wilcoxon rank sum test 

#data:  meta$invsimpson and meta$All 

# -3 INF_InjOTC -3 INF_NoTRMT -3 INF_OralOTC -3 NONINF_NoTRMT 0 INF_InjOTC 0 INF_NoTRMT
# -3 INF_NoTRMT    0.0796        -             -              -                -            -           
# -3 INF_OralOTC   0.7176        0.0620        -              -                -            -           
# -3 NONINF_NoTRMT 0.6702        0.1866        0.5532         -                -            -           
# 0 INF_InjOTC     0.1437        0.0065        0.1275         0.0620           -            -           
# 0 INF_NoTRMT     0.0273        0.0225        0.0190         0.0245           0.3929       -           
# 0 INF_OralOTC    0.1790        0.0318        0.0970         0.0943           1.0000       0.2677      
# 0 NONINF_NoTRMT  0.0428        0.0014        0.0264         0.0231           0.5556       0.7577      
# 7 INF_InjOTC     0.9382        0.0772        0.8121         0.5241           0.3128       0.0330      
# 7 INF_NoTRMT     0.3360        0.0356        0.2670         0.1458           0.6049       0.1613      
# 7 INF_OralOTC    0.2406        0.0122        0.1196         0.0642           0.9274       0.1331      
# 7 NONINF_NoTRMT  0.3196        0.0274        0.3285         0.1772           0.6970       0.0939      
# 0 INF_OralOTC 0 NONINF_NoTRMT 7 INF_InjOTC 7 INF_NoTRMT 7 INF_OralOTC
# -3 INF_NoTRMT    -             -               -            -            -            
# -3 INF_OralOTC   -             -               -            -            -            
# -3 NONINF_NoTRMT -             -               -            -            -            
# 0 INF_InjOTC     -             -               -            -            -            
# 0 INF_NoTRMT     -             -               -            -            -            
# 0 INF_OralOTC    -             -               -            -            -            
# 0 NONINF_NoTRMT  0.6993        -               -            -            -            
# 7 INF_InjOTC     0.1645        0.0947          -            -            -            
# 7 INF_NoTRMT     0.5135        0.2428          0.5288       -            -            
# 7 INF_OralOTC    0.6787        0.3154          0.2799       0.7959       -            
# 7 NONINF_NoTRMT  0.6216        0.3704          0.4082       0.8968       0.6965       

#P value adjustment method: none  

#Rename treatment groups in All and Treatment columns
meta$Treatment2 = meta$Treatment
meta$Treatment2 <- as.character(meta$Treatment2)
meta$Treatment2[meta$Treatment2 == 'INF_OralOTC'] <- "INF_OralOCT"
meta$Treatment2[meta$Treatment2 == 'INF_InjOTC'] <- "INF_InjOCT"
meta$All2 = meta$All
meta$All2 <- as.character(meta$All2)
meta$All2[meta$All2 == '-3 INF_OralOTC'] <- "-3 INF_OralOCT"
meta$All2[meta$All2 == '0 INF_OralOTC'] <- "0 INF_OralOCT"
meta$All2[meta$All2 == '7 INF_OralOTC'] <- "7 INF_OralOCT"
meta$All2[meta$All2 == '-3 INF_InjOTC'] <- "-3 INF_InjOCT"
meta$All2[meta$All2 == '0 INF_InjOTC'] <- "0 INF_InjOCT"
meta$All2[meta$All2 == '7 INF_InjOTC'] <- "7 INF_InjOCT"


#Generate a box and whisker plot of shannon (both shannon and inverse simpson diversity indices showed same trends; I chose to make a plot using shannon indices)
(shan <- ggplot(data = meta, aes(x=All2, y=shannon, group=All2, fill=Treatment2)) +
    geom_boxplot(position = position_dodge2(preserve = 'total')) +
    facet_wrap(~Treatment2, scales = 'free') +
    scale_y_continuous(name = "Shannon diversity") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          strip.text.x = element_text(size=14),
          axis.title.y = element_text(size=15)) +
    theme(legend.position = "none"))
# "free" within "facet_wrap" allows each plot to customize the scale to the specific data set (no forced scaling applied to all plots)
# "position = position_dodge2(preserve = 'total')" fixes the ggplot box width, making them wider, prevents narrow boxes from forming in the plot

#Save 'shan' as a .tiff for publication, 500dpi
ggsave("FS9_Shannon.tiff", plot=shan, width = 7, height = 7, dpi = 500, units =c("in"))

#Generate a box and whisker plot of inverse simpson 
(invsimp <- ggplot(data = meta, aes(x=All2, y=invsimpson, group=All2, fill=Treatment2)) +
    geom_boxplot(position = position_dodge2(preserve = 'total')) +
    facet_wrap(~Treatment2, scales = 'free') +
    scale_y_continuous(name = "Inverse Simpson diversity") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          strip.text.x = element_text(size=14),
          axis.title.y = element_text(size=15)) +
    theme(legend.position = "none"))

#Save 'invsimp' as a .tiff for publication, 500dpi
ggsave("FS9_InverseSimpson.tiff", plot=invsimp, width = 7, height = 7, dpi = 500, units =c("in"))
