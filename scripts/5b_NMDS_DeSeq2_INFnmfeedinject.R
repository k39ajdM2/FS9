#####################################################################################################
#FS9 NMDS plot and DeSeq2 (order level) plot combined into one figure for days 7, 11, 14 - Infected: nm vs feed or inject.
#By Mou, KT

#Purpose: Combine NMDS figure and PERMANOVA F magnitude of change figure for INFnm vs INFfeed vs INFinject on days 7, 11, 14 into one using cowplot package

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
library(DESeq2)

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

################################################ NMDS ################################################
#Generate phyloseq.FS9 object
meta <- read.csv("./data/FS9_metadata_INF_InjectFeedNM.csv", row.names = 1)
otu <- read.csv("./data/FS9.OTUtable.doubleton.csv", row.names=1)
dim(otu) #1405 168
head(otu[,165:168])

#Remove taxonomy from 'otu'
tax <- otu[,(167:168)] #Removed column 168 "Taxonomy" and 167 to include the row names
head(tax)
colnames(tax)[1] <- "delete" #Renamed column 1 (formerly 167) as "delete" which will later be deleted
head(tax)

#Modify 'otu' with only OTU count data
otu <- otu[,-168] #Remove column 168 "Taxonomy" to have only OTU data
head(otu[,165:167]) 
dim(otu) #1405 rows 167 columns

#Transpose 'otu' to match format of 'meta'
otu.trans <- t(otu) #Now rownames are sample names, columns are OTUs
head(otu.trans[,1:5])
class(meta) #dataframe
class(otu) #dataframe

#Merge 'otu' and 'meta' data frames
otu.meta <- merge(meta, otu.trans, by.x=0, by.y=0) 
#x=0 means match via rownames from 'meta'; y=0 means match via rownames from 'otu.trans'
head(otu.meta[,1:10])
dim(otu.meta) #72 1410 (doubletons removed)
rownames(otu.meta) <- otu.meta[,1] #Set column 1 as rownames for 'otu.meta'
class(otu.meta) #dataframe
otu.meta <- otu.meta[,-1]
rownames(otu.meta)

#Added an "All" column (combines 'Treatment' and 'Day' values) in new 'otu.meta2'
otu.meta2<- cbind(otu.meta) #Make second copy of otu.meta to use to include "All" column
colnames(otu.meta2)
otu.meta2$All <- with(otu.meta2, paste0(Day, sep="_", Treatment)) #Create "All" column with Day and Treatment combined
head(otu.meta2)
dim(otu.meta2) #72 1410 (doubletons removed)
head(otu.meta2[,1405:1410])
head(otu.meta2[,1:10])
otu.meta2<- otu.meta2[,c(1:4,1410,5:1409)] #Reorder columns to have "All" column after "Treatment" column
#write.csv(otu.meta2, file="FS9.otu.meta_all.doubleton.csv")

#Pull out metadata from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
meta2 <- otu.meta2[,c(1:5)] #Take columns 1-5 ("Day" to "All") from 'otu.meta2' to make 'meta2'
head(meta2)
dim(meta2) #72 5 (doubletons removed)

#Create SAM metadata table phyloseq object
SAM = sample_data(meta2, errorIfNULL = TRUE)
head(SAM)
dim(SAM) #72 5

#Pull out OTU data from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
tail(otu.meta2[,1405:1410])
otu2 <- otu.meta2[,c(6:1410)] #Select OTU columns to create 'otu2' dataframe
head(otu2[,1:10])
dim(otu2) #72 1405 (doubletons removed)
otu2.trans <- t(otu2) #Transpose otu.all to have OTUs as rownames, sample names as column names
head(otu2.trans[,1:10])
dim(otu2.trans) #1405 72 (doubletons removed)

#Merge 'tax' back into 'otu2.trans' for correct format and taxons
head(tax)
otu2.tax <- merge(otu2.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu2.tax) #1405 75 (doubletons removed)
head(otu2.tax[,1:10])
head(otu2.tax[,70:75])
row.names(otu2.tax) <- otu2.tax[,1] #Set first row of 'otu2.tax' as rownames
head(otu2.tax[,1:5])
otu2.tax <- otu2.tax[,-1] #Remove first row aka extraneous "Row.names" column from 'otu2.tax'
head(otu2.tax[,1:5])

#Split again
dim(otu2.tax) #1405 74 (doubletons removed)
head(otu2.tax[,70:74])
otu2.notax <- otu2.tax[,1:72] #take rows 1-72 to make new dataframe 'otu2.notax' (73 is "delete" column, 74 is "Taxonomy" column)
dim(otu2.notax) #1405 72 (doubletons removed)
head(otu2.notax[,1:5])
head(otu2.notax[,70:72])
class(otu2.notax) #dataframe
otu2.notax <- as.matrix(otu2.notax) #turn 'otu2.notax' into a matrix class
class(otu2.notax) #matrix array
otu2.notax
otu2.notax.trans <- t(otu2.notax)
head(otu2.notax.trans[,1:10])

#Create OTU table phyloseq object
OTU = otu_table(otu2.notax.trans, taxa_are_rows = FALSE)
head(OTU[,1:10])
dim(OTU) #72 1405 (doubletons removed)
class(OTU)
OTU2 <- prune_taxa(taxa_sums(OTU) > 0, OTU)
class(OTU2)
head(OTU2[,1:10])
taxa_sums(OTU2)
dim(OTU2) #72 1042 (doubletons removed)

#Edit taxonomy
dim(otu2.tax) #1405 72 (doubletons removed)
head(otu2.tax[,70:72])
tax2 <- separate(data = otu2.tax, 
                 col = Taxonomy, 
                 into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#"separate" function separates "Taxonomy" column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax2) #notice that Species column is blank
dim(tax2) #1405 80 (doubletons removed)
head(tax2[,70:80])
tax2.kg <- tax2[,74:79] #Keep only taxonomy columns "Kingdom" to "Genus"
head(tax2.kg)
dim(tax2.kg) #1405 6 (doubletons removed)
class(tax2.kg) #dataframe
tax2.kg <- as.matrix(tax2.kg)
class(tax2.kg) #matrix array
tax2.kg

#Create TAX taxonomy table phyloseq object
TAX = tax_table(tax2.kg)
head(TAX)
dim(TAX) #1405 6 (doubletons removed)

#Create phyloseq object containing taxonomy, metadata, and otu table
phyloseq.FS9 <- phyloseq(OTU2, SAM, TAX) #prune out any OTUs that have total to 0 in all samples
phyloseq.FS9
#otu_table()   OTU Table:         [ 1042 taxa and 72 samples ]
#sample_data() Sample Data:       [ 72 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 1042 taxa by 6 taxonomic ranks ]

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
  labs(color="Treatment")
nmdsplot_treatment2

##################################################### DESeq2 ######################################################################

otu <- import_mothur(mothur_shared_file = './data/stability.outdoubletons.abund.opti_mcc.shared') #use unrarified data
taxo <- import_mothur(mothur_constaxonomy_file = './data/stability.outdoubletons.abund.opti_mcc.0.03.cons.taxonomy')
meta <- read.table(file = './data/FS9_metadata_INF_InjectFeedNM.csv', sep = ',', header = TRUE)

#Organize meta file
rownames(meta) <- meta$Sample
meta <- meta[,-1] #remove Sample column
meta$Set <- paste(meta$Day, meta$Treatment, sep = '_') #Make a set column

#Make phyloseq object SRD129 (combine taxonomy, OTU, and metadata)
phy_meta <- sample_data(meta) 
FS9 <- phyloseq(otu, taxo)
FS9 <- merge_phyloseq(FS9, phy_meta)   # combines the metadata with this phyloseq object
colnames(tax_table(FS9))
colnames(tax_table(FS9)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')

FS9
#phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 2039 taxa and 73 samples ]
#sample_data() Sample Data:       [ 73 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 2039 taxa by 6 taxonomic ranks ]

#Prune
FS9 <- prune_samples(sample_sums(FS9) > 1300, FS9)  # This removes samples that have fewer than 1300 sequences associated with them.
#Jules used 1000. You used 1300 because Day 0 Treatment group 2 (room 2) only has 5 samples and subsampling at 2000 would take out 
#one of the samples.
FS9 <- prune_taxa(taxa_sums(FS9) > 10, FS9)        # removes OTUs that occur less than 10 times globally
tax_table(FS9) [1:5, 1:6] #see what's in tax_table first 5 rows, first 6 columns

# If you want to group OTUs uncomment the tax_glom() line and select your desired taxrank
# right now all these analysis are done at the OTU level.

#Grouping OTUs by desired taxonomic level using the tax_glom function 
FS9.order <- tax_glom(FS9, taxrank = "Order")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 

##################################################### Day 14 ######################################################################

sample_data(FS9.order)

FS9.D14 <- subset_samples(FS9.order, Day == 'D14')
sample_sums(FS9.D14)
colnames(otu_table(FS9.D14)) #check on all the sample names
FS9.D14 <- prune_taxa(taxa_sums(FS9.D14) > 1, FS9.D14)
#if taxa_sums is >1, then it will print that out in FS9.D14 object and not include anything with <1.
rowSums(FS9.D14@otu_table)

######### 1. Day 14 INFinject vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D14_INFnm")
#14
sum(meta$Set == "D14_INFinject")
#19

sample_data(FS9.D14)$Set <- factor(sample_data(FS9.D14)$Set,
                                   levels =c("D14_INFinject", 'D14_INFnm', "D14_INFfeed"))
FS9.D14.De <- phyloseq_to_deseq2(FS9.D14, ~ Set)
FS9.D14.De <- DESeq(FS9.D14.De, test = "Wald", fitType = "parametric")
FS9.D14.De$Set
resultsNames(FS9.D14.De)
#[1] "Intercept"        "Set_D14_INFnm_vs_D14_INFinject"   "Set_D14_INFfeed_vs_D14_INFinject"
res.D14.ji = lfcShrink(FS9.D14.De, coef = "Set_D14_INFnm_vs_D14_INFinject", type = 'apeglm')
sigtab.D14.ji = res.D14.ji[which(res.D14.ji$padj < .05), ]
sigtab.D14.ji = cbind(as(sigtab.D14.ji, "data.frame"), as(tax_table(FS9.D14)[rownames(sigtab.D14.ji), ], "matrix"))
sigtab.D14.ji$newp <- format(round(sigtab.D14.ji$padj, digits = 3), scientific = TRUE)
sigtab.D14.ji$Treatment <- ifelse(sigtab.D14.ji$log2FoldChange >=0, "Depleted in INFinject", "Enriched in INFinject")
head(sigtab.D14.ji) 

deseq.D14.ji <- 
  ggplot(sigtab.D14.ji, aes(x=reorder(rownames(sigtab.D14.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.ji), y=3, label = paste(Order)), size=5, fontface = 'italic')+ labs(x="Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFinject Group Relative to INFnm\n in Fecal Microbiota on Day 14')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c('Enriched INFinject'='#00BA38',
                               'Depleted in INFinject'='#F8766D'))
deseq.D14.ji

#Add OTU and comparisons columns
sigtab.D14.ji
sigtab.D14.ji$OTU <- rownames(sigtab.D14.ji)
sigtab.D14.ji
sigtab.D14.ji$comp <- 'D14_INFnm_vs_INFinject'

#Create final significant comparisons table
final.sigtab.ji <- sigtab.D14.ji


######### 2. Day 14 INFfeed vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D14_INFnm")
#14
sum(meta$Set == "D14_INFfeed")
#13

#Extract results from a DESeq analysis, organize table
sample_data(FS9.D14)$Set <- factor(sample_data(FS9.D14)$Set,
                                   levels =c( "D14_INFfeed", 'D14_INFnm',"D14_INFinject"))
FS9.D14.De <- phyloseq_to_deseq2(FS9.D14, ~ Set)
FS9.D14.De <- DESeq(FS9.D14.De, test = "Wald", fitType = "parametric")
FS9.D14.De$Set
resultsNames(FS9.D14.De)
#"Set_D14_INFnm_vs_D14_INFfeed"     "Set_D14_INFinject_vs_D14_INFfeed"
res.D14.oi = lfcShrink(FS9.D14.De, coef = "Set_D14_INFnm_vs_D14_INFfeed", type = 'apeglm')
sigtab.D14.oi = res.D14.oi[which(res.D14.oi$padj < .05), ]
sigtab.D14.oi = cbind(as(sigtab.D14.oi, "data.frame"), as(tax_table(FS9.D14)[rownames(sigtab.D14.oi), ], "matrix"))
format(sigtab.D14.oi$padj, scientific = TRUE)
sigtab.D14.oi$newp <- format(round(sigtab.D14.oi$padj, digits = 3), scientific = TRUE)
sigtab.D14.oi$Treatment <- ifelse(sigtab.D14.oi$log2FoldChange >=0, "Depleted in INFfeed", "Enriched in INFfeed")
head(sigtab.D14.oi) 

deseq.D14.oi <- 
  ggplot(sigtab.D14.oi, aes(x=reorder(rownames(sigtab.D14.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.oi), y=3, label = paste(Order)), size=5, fontface = 'italic')+ labs(x="Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFfeed Group Relative to INFnm\n in Fecal Microbiota on Day 14')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c('Enriched INFfeed'='#619CFF',
                               'Depleted in INFfeed'='#F8766D'))
deseq.D14.oi

#Add OTU and comparisons columns
sigtab.D14.oi
sigtab.D14.oi$OTU <- rownames(sigtab.D14.oi)
sigtab.D14.oi
sigtab.D14.oi$comp <- 'D14_INFfeed_vs_INFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab.ji, sigtab.D14.oi)

#Figure 4 original plot containing all significant orders with log2-fold changes > 0.25 for Q2 DeSeq2 results
deseqfinalplot <- ggplot(final.sigtab, aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity', position = position_dodge(preserve = 'single')) + geom_text(y=2, label = paste(final.sigtab$Order), size=5, fontface = 'italic')+ labs(x="Order")+
  coord_flip() +
  scale_fill_manual(values = c('Depleted in INFfeed'='lightskyblue', 'Depleted in INFinject'='lightgreen'), name= "Treatment relative to INFnm") +
  theme_bw() +
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12), 
        axis.title.x=element_text(size = 12),
        axis.title.y =element_text(size = 12),
        legend.text = element_text(size=11), 
        legend.title = element_text(size=12))
deseqfinalplot
#light green for depleted in INFinject
#light blue for INFfeed
#color name options for ggplot: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf

#Cowplot of NMDS with DESeq2
Fig4ab<- plot_grid(nmdsplot_treatment2, deseqfinalplot, labels = c('A', 'B'), label_size = 18, ncol=1)
Fig4ab

ggsave("INFnm_inject_feed_NMDS_DESeq2.tiff", plot=Fig4ab, width = 10, height = 9, dpi = 500, units =c("in"))
