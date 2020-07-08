#####################################################################################################
#FS9 Nasal DESeq2
#Kathy Mou

#Purpose: This code uses DESeq2 package to identify nasal microbial genera that were differentially 
#abundant between the three groups on the three days sampled

#Files needed:


#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory (either on local or network drive, Mac or PC)
#Mac
setwd("~/Desktop/FS9/FS9_RWorkspace")

#Load library packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("XML")
BiocManager::install("annotate")
BiocManager::install("genefilter")
BiocManager::install("geneplotter")
BiocManager::install("DESeq2")
#Couldn't get DESeq2 library to install

#Load library packages
library(DESeq2)
library(phyloseq)
library(ggplot2)
library("wesanderson")
library(plotly)
library(gapminder)
library("ggsci")
library("apeglm")

sessionInfo()

# Jules note:  I reccommend against doing things like this
# best to load the data from csv files etc
# for example, I have no idea what commands generated the .RData files you
# are loading here so when you're having issues I need to start from the 
# raw data files


#Load image file 
load("FS9_DESeq2.RData")




#Save image file
save.image(file="FS9_DESeq2.RData")

#Additional notes
#Jules says plots describes the log-fold changes seen in differential abundance plots as
#enriched for "x" taxa in "x group"
#You can add genus to your plots, and use those as supplemental figures


#Set plots to have gray background with white gridlines
theme_set(theme_gray())

########################################################################################################

####### PREPARING OBJECTS FOR DESEQ2 ANALYSIS ########
# Jules note: These lines of code don't work, it's best to make it
# so your scripts will run line by line from beginning to end without 
# needing to jump around.

# I added the paths needed but even this is not very optimal.
# I like to use Rprojects, that way your paths are all relative to 
# the project and as long as you move the whole project directory it 
# will work in and location 

#Load files
otu <- import_mothur(mothur_shared_file = '/ssd/kathy.mou/FS9_RWorkspace/stability.outsingletons.abund.opti_mcc.shared') #use unrarified data
taxo <- import_mothur(mothur_constaxonomy_file = '/ssd/kathy.mou/FS9_RWorkspace/stability.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
meta <- read.table(file = '/ssd/kathy.mou/FS9_RWorkspace/FS9_metadata.csv', sep = ',', header = TRUE)

#Organize meta file
rownames(meta) <- meta$Sample
meta <- meta[,-1] #remove Sample column
meta <- meta[,-5] #remove All column and replace with new Set column (below)
meta$Set <- paste(meta$Day, meta$Treatment, sep = '_')

#Make phyloseq object SRD129 (combine taxonomy, OTU, and metadata)
phy_meta <- sample_data(meta) 
FS9 <- phyloseq(otu, taxo)
FS9 <- merge_phyloseq(FS9, phy_meta)   # combines the metadata with this phyloseq object
colnames(tax_table(FS9))
colnames(tax_table(FS9)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')

FS9

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
FS9.phylum <- tax_glom(FS9, taxrank = "Phylum")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 

######################################## PRIMARY COMPARISONS TO MAKE ############################################################

# NMDS plot showed that disperion is different between days, so I subsetted by day

# Important comparisons to make (significant changes in beta diversity between treatments): 
# Day -3 
# INF_NoTRMT vs INF_InjOTC 
# INF_NoTRMT vs INF_OralOTC
# NONINF_NoTRMT vs INF_NoTRMT 
# NONINF_NoTRMT vs INF_OralOTC

# Day 7
# INF_NoTRMT vs INF_InjOTC
# INF_NoTRMT vs INF_OralOTC
# INF_InjOTC vs INF_OralOTC
# INF_InjOTC vs NONINF_NoTRMT


# Other comparisons to make (no significant changes in beta diversity between treatments): 
# Day -3
# INF_OralOTC vs INF_InjOTC
# NONINF_NoTRMT vs INF_InjOTC

# Day 0
# INF_InjOTC vs INF_NoTRMT
# INF_OralOTC vs INF_InjOTC
# INF_OralOTC vs INF_NoTRMT
# NONINF_NoTRMT vs INF_InjOTC
# NONINF_NoTRMT vs INF_NoTRMT
# NONINF_NoTRMT vs INF_OralOTC

# Day 7
# NONINF_NoTRMT vs INF_NoTRMT
# NONINF_NoTRMT vs INF_OralOTC

##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (ORDER LEVEL) BETWEEN TREATMENTS ###################################

######################################################### Day -3 #########################################################

sample_data(FS9.order)

FS9.DNEG3 <- subset_samples(FS9.order, Day == '-3')
sample_sums(FS9.DNEG3)
colnames(otu_table(FS9.DNEG3)) #check on all the sample names
FS9.DNEG3 <- prune_taxa(taxa_sums(FS9.DNEG3) > 1, FS9.DNEG3)
#if taxa_sums is >1, then it will print that out in FS9.DNEG3 object and not include anything with <1.
rowSums(FS9.DNEG3@otu_table)

#Look at what Set is
sample_data(FS9.DNEG3)

FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing

######### Day -3 INF_NoTRMT vs INF_InjOTC ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC= J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "-3_INF_NoTRMT")
#20
sum(meta$Set == "-3_INF_InjOTC")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set

# Jules Add
resultsNames(FS9.DNEG3.De)

res.DNEG3.ji = results(FS9.DNEG3.De, contrast=c("Set", "-3_INF_InjOTC", "-3_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')

res.DNEG3.ji = lfcShrink(FS9.DNEG3.De, res=res.DNEG3.ji, coef = "2", type = 'apeglm')
#Error in lfcShrink(FS9.DNEG3.De, res = res.DNEG3.ji, coef = 2, type = "apeglm") : 
#'coef' should specify same coefficient as in results 'res'
#I tried reordering the contrasts so that it was listed as:
#contrast=c("Set", "-3_INF_NoTRMT", "-3_INF_InjOTC") to match vector 2 (from resultsNames(FS9.DNEG3.De))
#but I still got the same error message

res.DNEG3.ji = lfcShrink(FS9.DNEG3.De, res=res.DNEG3.ji, coef = "Set_.3_INF_NoTRMT_vs_.3_INF_InjOTC", type = 'apeglm')
#Error in lfcShrink(FS9.DNEG3.De, res = res.DNEG3.ji, coef = "Set_.3_INF_NoTRMT_vs_.3_INF_InjOTC",  : 
#'coef' should specify same coefficient as in results 'res'
#Same error message

sigtab.DNEG3.ji = res.DNEG3.ji[which(res.DNEG3.ji$padj < .05), ]
sigtab.DNEG3.ji = cbind(as(sigtab.DNEG3.ji, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.ji), ], "matrix"))
format(sigtab.DNEG3.ji$padj, scientific = TRUE)
sigtab.DNEG3.ji$newp <- format(round(sigtab.DNEG3.ji$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.ji$Treatment <- ifelse(sigtab.DNEG3.ji$log2FoldChange >=0, "INF_InjOTC", "INF_NoTRMT")

#Summarize sigtab.DNEG3.ji
sum.sigtab.DNEG3.ji <- summary(sigtab.DNEG3.ji)
sum.sigtab.DNEG3.ji



# your scale_fill_manual is the problem
# This is the same code block as yours except I commented out the scale fill manual

deseq.DNEG3.ji_JULES <- 
  ggplot(sigtab.DNEG3.ji, aes(x=reorder(rownames(sigtab.DNEG3.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.ji), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_InjOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) #+
  # scale_fill_manual(labels = c("INF_InjOTC", "INF_NoTRMT"), values = c('#E69F00', '#CC0066'))
deseq.DNEG3.ji_JULES

# you were expecting the labels to match the levels in your Treatment factor, 
# but when you pass in a 'labels' vector is overrides the original labels.
# if you want to do the matching thing you were expecting you need to pass in
# a named vector, like I do in the next two blocks

deseq.DNEG3.ji_JULES2 <- 
  ggplot(sigtab.DNEG3.ji, aes(x=reorder(rownames(sigtab.DNEG3.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.ji), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_InjOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
scale_fill_manual(values = c(INF_InjOTC='#E69F00', INF_NoTRMT='#CC0066'))

deseq.DNEG3.ji_JULES2


deseq.DNEG3.ji_JULES3 <- 
  ggplot(sigtab.DNEG3.ji, aes(x=reorder(rownames(sigtab.DNEG3.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.ji), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_InjOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INF_NoTRMT='#E69F00', INF_InjOTC='#CC0066'))

deseq.DNEG3.ji_JULES3


######### END JULES #######

#ggplot
deseq.DNEG3.ji <- 
  ggplot(sigtab.DNEG3.ji, aes(x=reorder(rownames(sigtab.DNEG3.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.ji), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_InjOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INF_NoTRMT='#E69F00', INF_InjOTC='#CC0066'))
deseq.DNEG3.ji

#Add OTU and comparisons columns
sigtab.DNEG3.ji
sigtab.DNEG3.ji$OTU <- rownames(sigtab.DNEG3.ji)
sigtab.DNEG3.ji
sigtab.DNEG3.ji$comp <- 'DNEG3_INF_InjOTCvsINF_NoTRMT'

#Create final significant comparisons table
final.sigtab <- sigtab.DNEG3.ji



######### Day -3 INF_NoTRMT vs INF_OralOTC ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "-3_INF_NoTRMT")
#20
sum(meta$Set == "-3_INF_OralOTC")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
res.DNEG3.oi = results(FS9.DNEG3.De, contrast=c("Set", "-3_INF_OralOTC", "-3_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.DNEG3.De, contrast=c("Set","-3_INF_OralOTC", "-3_INF_NoTRMT")) 
sigtab.DNEG3.oi = res.DNEG3.oi[which(res.DNEG3.oi$padj < .05), ]
sigtab.DNEG3.oi = cbind(as(sigtab.DNEG3.oi, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.oi), ], "matrix"))
format(sigtab.DNEG3.oi$padj, scientific = TRUE)
sigtab.DNEG3.oi$newp <- format(round(sigtab.DNEG3.oi$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.oi$Treatment <- ifelse(sigtab.DNEG3.oi$log2FoldChange >=0, "INF_OralOTC", "INF_NoTRMT")

#Summarize sigtab.DNEG3.oi
sum.sigtab.DNEG3.oi <- summary(sigtab.DNEG3.oi)
sum.sigtab.DNEG3.oi

#ggplot
deseq.DNEG3.oi <- ggplot(sigtab.DNEG3.oi, aes(x=reorder(rownames(sigtab.DNEG3.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.oi), y=-1, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_OralOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_OralOTC", "INF_NoTRMT"), values = c('#999999', '#CC0066'))
  scale_fill_manual(values = c(INF_OralOTC='#999999', INF_NoTRMT='#CC0066'))
deseq.DNEG3.oi

#Add OTU and comparisons columns
sigtab.DNEG3.oi
sigtab.DNEG3.oi$OTU <- rownames(sigtab.DNEG3.oi)
sigtab.DNEG3.oi
sigtab.DNEG3.oi$comp <- 'DNEG3_INF_OralOTCvsINF_NoTRMT'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.DNEG3.oi, final.sigtab)



######### Day -3 NONINF_NoTRMT vs INF_NoTRMT  ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "-3_NONINF_NoTRMT")
#20
sum(meta$Set == "-3_INF_NoTRMT")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
res.DNEG3.in = results(FS9.DNEG3.De, contrast=c("Set", "-3_INF_NoTRMT", "-3_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.DNEG3.De, contrast=c("Set","-3_INF_NoTRMT", "-3_NONINF_NoTRMT")) 
sigtab.DNEG3.in = res.DNEG3.in[which(res.DNEG3.in$padj < .05), ]
sigtab.DNEG3.in = cbind(as(sigtab.DNEG3.in, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.in), ], "matrix"))
format(sigtab.DNEG3.in$padj, scientific = TRUE)
sigtab.DNEG3.in$newp <- format(round(sigtab.DNEG3.in$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.in$Treatment <- ifelse(sigtab.DNEG3.in$log2FoldChange >=0, "INF_NoTRMT", "NONINF_NoTRMT")

#Summarize sigtab.DNEG3.in
sum.sigtab.DNEG3.in <- summary(sigtab.DNEG3.in)
sum.sigtab.DNEG3.in

#ggplot
deseq.DNEG3.in <- ggplot(sigtab.DNEG3.in, aes(x=reorder(rownames(sigtab.DNEG3.in), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.in), y=2, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_NoTRMT Group Relative to NONINF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_NoTRMT", "NONINF_NoTRMT"), values = c('#CC0066', '#56B4E9'))
  scale_fill_manual(values = c(INF_NoTRMT='#CC0066', NONINF_NoTRMT='#56B4E9'))
deseq.DNEG3.in

#Add OTU and comparisons columns
sigtab.DNEG3.in
sigtab.DNEG3.in$OTU <- rownames(sigtab.DNEG3.in)
sigtab.DNEG3.in
sigtab.DNEG3.in$comp <- 'DNEG3_INF_NoTRMTvsNONINF_NoTRMT'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.DNEG3.in)



######### Day -3 NONINF_NoTRMT vs INF_OralOTC  ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "-3_NONINF_NoTRMT")
#20
sum(meta$Set == "-3_INF_OralOTC")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
res.DNEG3.on = results(FS9.DNEG3.De, contrast=c("Set", "-3_INF_OralOTC", "-3_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.DNEG3.De, contrast=c("Set","-3_INF_OralOTC", "-3_NONINF_NoTRMT")) 
sigtab.DNEG3.on = res.DNEG3.on[which(res.DNEG3.on$padj < .05), ]
sigtab.DNEG3.on = cbind(as(sigtab.DNEG3.on, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.on), ], "matrix"))
format(sigtab.DNEG3.on$padj, scientific = TRUE)
sigtab.DNEG3.on$newp <- format(round(sigtab.DNEG3.on$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.on$Treatment <- ifelse(sigtab.DNEG3.on$log2FoldChange >=0, "INF_OralOTC", "NONINF_NoTRMT")

#Summarize sigtab.DNEG3.on
sum.sigtab.DNEG3.on <- summary(sigtab.DNEG3.on)
sum.sigtab.DNEG3.on

#ggplot
deseq.DNEG3.on <- ggplot(sigtab.DNEG3.on, aes(x=reorder(rownames(sigtab.DNEG3.on), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.on), y=-1, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_OralOTC Group Relative to NONINF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_OralOTC", "NONINF_NoTRMT"), values = c('#999999', '#56B4E9'))
  scale_fill_manual(values = c(INF_OralOTC='#999999', NONINF_NoTRMT='#56B4E9'))
deseq.DNEG3.on

#Add OTU and comparisons columns
sigtab.DNEG3.on
sigtab.DNEG3.on$OTU <- rownames(sigtab.DNEG3.on)
sigtab.DNEG3.on
sigtab.DNEG3.on$comp <- 'DNEG3_INF_OralOTCvsNONINF_NoTRMT'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.DNEG3.on)


##################################################### Day 7 ######################################################################

sample_data(FS9.order)

FS9.D7 <- subset_samples(FS9.order, Day == '7')
sample_sums(FS9.D7)
colnames(otu_table(FS9.D7)) #check on all the sample names
FS9.D7 <- prune_taxa(taxa_sums(FS9.D7) > 1, FS9.D7)
#if taxa_sums is >1, then it will print that out in FS9.D7 object and not include anything with <1.
rowSums(FS9.D7@otu_table)

#Look at what Set is
sample_data(FS9.D7)

FS9.D7.De <- phyloseq_to_deseq2(FS9.D7, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.D7.De <- DESeq(FS9.D7.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing

######### Day 7 INF_InjOTC vs INF_NoTRMT  ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC= J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "7_INF_InjOTC")
#26
sum(meta$Set == "7_INF_NoTRMT")
#20

#Extract results from a DESeq analysis, organize table
FS9.D7.De$Set
res.D7.ji = results(FS9.D7.De, contrast=c("Set", "7_INF_InjOTC", "7_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.D7.De, contrast=c("Set","7_INF_InjOTC", "7_INF_NoTRMT")) 
sigtab.D7.ji = res.D7.ji[which(res.D7.ji$padj < .05), ]
sigtab.D7.ji = cbind(as(sigtab.D7.ji, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7.ji), ], "matrix"))
format(sigtab.D7.ji$padj, scientific = TRUE)
sigtab.D7.ji$newp <- format(round(sigtab.D7.ji$padj, digits = 3), scientific = TRUE)
sigtab.D7.ji$Treatment <- ifelse(sigtab.D7.ji$log2FoldChange >=0, "INF_InjOTC", "INF_NoTRMT")

#Summarize sigtab.D7.ji
sum.sigtab.D7.ji <- summary(sigtab.D7.ji)
sum.sigtab.D7.ji

#ggplot
deseq.D7.ji <- ggplot(sigtab.D7.ji, aes(x=reorder(rownames(sigtab.D7.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.ji), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_InjOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_InjOTC", "INF_NoTRMT"), values = c('#E69F00', '#CC0066'))
  scale_fill_manual(values = c(INF_InjOTC='#E69F00', INF_NoTRMT='#CC0066'))
deseq.D7.ji

#Add OTU and comparisons columns
sigtab.D7.ji
sigtab.D7.ji$OTU <- rownames(sigtab.D7.ji)
sigtab.D7.ji
sigtab.D7.ji$comp <- 'D7_INF_InjOTCvsINF_NoTRMT'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.D7.ji, final.sigtab)



######### Day 7 INF_OralOTC vs INF_NoTRMT ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "7_INF_NoTRMT")
#20
sum(meta$Set == "7_INF_OralOTC")
#17

#Extract results from a DESeq analysis, organize table
FS9.D7.De$Set
res.D7.oi = results(FS9.D7.De, contrast=c("Set", "7_INF_OralOTC", "7_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.D7.De, contrast=c("Set","7_INF_OralOTC", "7_INF_NoTRMT")) 
sigtab.D7.oi = res.D7.oi[which(res.D7.oi$padj < .05), ]
sigtab.D7.oi = cbind(as(sigtab.D7.oi, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7.oi), ], "matrix"))
format(sigtab.D7.oi$padj, scientific = TRUE)
sigtab.D7.oi$newp <- format(round(sigtab.D7.oi$padj, digits = 3), scientific = TRUE)
sigtab.D7.oi$Treatment <- ifelse(sigtab.D7.oi$log2FoldChange >=0, "INF_OralOTC", "INF_NoTRMT")

#Summarize sigtab.D7.oi
sum.sigtab.D7.oi <- summary(sigtab.D7.oi)
sum.sigtab.D7.oi

#ggplot
deseq.D7.oi <- ggplot(sigtab.D7.oi, aes(x=reorder(rownames(sigtab.D7.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.oi), y=-3, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_OralOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_OralOTC", "INF_NoTRMT"), values = c('#999999', '#CC0066'))
  scale_fill_manual(values = c(INF_OralOTC='#999999', INF_NoTRMT='#CC0066'))
deseq.D7.oi

#Add OTU and comparisons columns
sigtab.D7.oi
sigtab.D7.oi$OTU <- rownames(sigtab.D7.oi)
sigtab.D7.oi
sigtab.D7.oi$comp <- 'D7_INF_OralOTCvsINF_NoTRMT'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D7.oi)




######### Day 7 INF_OralOTC vs INF_InjOTC ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "7_INF_InjOTC")
#26
sum(meta$Set == "7_INF_OralOTC")
#17

#Extract results from a DESeq analysis, organize table
FS9.D7.De$Set
res.D7.oj = results(FS9.D7.De, contrast=c("Set", "7_INF_OralOTC", "7_INF_InjOTC"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.D7.De, contrast=c("Set","7_INF_OralOTC", "7_INF_InjOTC")) 
sigtab.D7.oj = res.D7.oj[which(res.D7.oj$padj < .05), ]
sigtab.D7.oj = cbind(as(sigtab.D7.oj, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7.oj), ], "matrix"))
format(sigtab.D7.oj$padj, scientific = TRUE)
sigtab.D7.oj$newp <- format(round(sigtab.D7.oj$padj, digits = 3), scientific = TRUE)
sigtab.D7.oj$Treatment <- ifelse(sigtab.D7.oj$log2FoldChange >=0, "INF_OralOTC", "INF_InjOTC")

#Summarize sigtab.D7.oj
sum.sigtab.D7.oj <- summary(sigtab.D7.oj)
sum.sigtab.D7.oj

#ggplot
deseq.D7.oj <- ggplot(sigtab.D7.oj, aes(x=reorder(rownames(sigtab.D7.oj), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.oj), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_OralOTC Group Relative to INF_InjOTC in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_OralOTC", "INF_InjOTC"), values = c('#999999', '#E69F00'))
  scale_fill_manual(values = c(INF_OralOTC='#999999', INF_InjOTC='#E69F00'))
deseq.D7.oj

#Add OTU and comparisons columns
sigtab.D7.oj
sigtab.D7.oj$OTU <- rownames(sigtab.D7.oj)
sigtab.D7.oj
sigtab.D7.oj$comp <- 'D7_INF_OralOTCvsINF_InjOTC'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D7.oj)



######### Day 7 INF_InjOTC vs NONINF_NoTRMT  ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC= J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "7_INF_InjOTC")
#26
sum(meta$Set == "7_NONINF_NoTRMT")
#17

#Extract results from a DESeq analysis, organize table
FS9.D7.De$Set
res.D7.jn = results(FS9.D7.De, contrast=c("Set", "7_INF_InjOTC", "7_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.D7.De, contrast=c("Set","7_INF_InjOTC", "7_NONINF_NoTRMT")) 
sigtab.D7.jn = res.D7.jn[which(res.D7.jn$padj < .05), ]
sigtab.D7.jn = cbind(as(sigtab.D7.jn, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7.jn), ], "matrix"))
format(sigtab.D7.jn$padj, scientific = TRUE)
sigtab.D7.jn$newp <- format(round(sigtab.D7.jn$padj, digits = 3), scientific = TRUE)
sigtab.D7.jn$Treatment <- ifelse(sigtab.D7.jn$log2FoldChange >=0, "INF_InjOTC", "NONINF_NoTRMT")

#Summarize sigtab.D7.jn
sum.sigtab.D7.jn <- summary(sigtab.D7.jn)
sum.sigtab.D7.jn

#ggplot
deseq.D7.jn <- ggplot(sigtab.D7.jn, aes(x=reorder(rownames(sigtab.D7.jn), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.jn), y=1, label = paste(Phylum,Order, sep = ' ')), size=5)+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INF_InjOTC Group Relative to NONINF_NoTRMT in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_InjOTC", "NONINF_NoTRMT"), values = c('#E69F00', '#56B4E9'))
  scale_fill_manual(values = c(INF_InjOTC='#E69F00', NONINF_NoTRMT='#56B4E9'))
deseq.D7.jn

#Add OTU and comparisons columns
sigtab.D7.jn
sigtab.D7.jn$OTU <- rownames(sigtab.D7.jn)
sigtab.D7.jn
sigtab.D7.jn$comp <- 'D7_INF_InjOTCvsNONINF_NoTRMT'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.D7.jn, final.sigtab)



#######################################################################################################

#write csv
write.csv(final.sigtab, file= "FS9_FinalDiffAbund_SignificantDays_Order.csv")

#######################################################################################################






##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (PHYLUM LEVEL) BETWEEN TREATMENTS ###################################

######################################################### Day -3 #########################################################

sample_data(FS9.phylum)

FS9.DNEG3.p <- subset_samples(FS9.phylum, Day == '-3')
sample_sums(FS9.DNEG3.p)
colnames(otu_table(FS9.DNEG3.p)) #check on all the sample names
FS9.DNEG3.p <- prune_taxa(taxa_sums(FS9.DNEG3.p) > 1, FS9.DNEG3.p)
#if taxa_sums is >1, then it will print that out in FS9.DNEG3.p object and not include anything with <1.
rowSums(FS9.DNEG3.p@otu_table)

#Look at what Set is
sample_data(FS9.DNEG3.p)

FS9.DNEG3.p.De <- phyloseq_to_deseq2(FS9.DNEG3.p, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.DNEG3.p.De <- DESeq(FS9.DNEG3.p.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing

######### Day -3 INF_NoTRMT vs INF_InjOTC ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC= J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "-3_INF_NoTRMT")
#20
sum(meta$Set == "-3_INF_InjOTC")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.p.De$Set
res.DNEG3.p.ji = results(FS9.DNEG3.p.De, contrast=c("Set", "-3_INF_InjOTC", "-3_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.DNEG3.p.De, contrast=c("Set","-3_INF_InjOTC", "-3_INF_NoTRMT")) 
sigtab.DNEG3.p.ji = res.DNEG3.p.ji[which(res.DNEG3.p.ji$padj < .05), ]
sigtab.DNEG3.p.ji = cbind(as(sigtab.DNEG3.p.ji, "data.frame"), as(tax_table(FS9.DNEG3.p)[rownames(sigtab.DNEG3.p.ji), ], "matrix"))
format(sigtab.DNEG3.p.ji$padj, scientific = TRUE)
sigtab.DNEG3.p.ji$newp <- format(round(sigtab.DNEG3.p.ji$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.p.ji$Treatment <- ifelse(sigtab.DNEG3.p.ji$log2FoldChange >=0, "INF_InjOTC", "INF_NoTRMT")

#Summarize sigtab.DNEG3.p.ji
sum.sigtab.DNEG3.p.ji <- summary(sigtab.DNEG3.p.ji)
sum.sigtab.DNEG3.p.ji

#ggplot
deseq.DNEG3.p.ji <- ggplot(sigtab.DNEG3.p.ji, aes(x=reorder(rownames(sigtab.DNEG3.p.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p.ji), y=-2, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phyla in INF_InjOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_InjOTC", "INF_NoTRMT"), values = c('#E69F00', '#CC0066'))
  scale_fill_manual(values = c(INF_InjOTC='#E69F00', INF_NoTRMT='#CC0066'))
deseq.DNEG3.p.ji

#Add OTU and comparisons columns
sigtab.DNEG3.p.ji
sigtab.DNEG3.p.ji$OTU <- rownames(sigtab.DNEG3.p.ji)
sigtab.DNEG3.p.ji
sigtab.DNEG3.p.ji$comp <- 'DNEG3_INF_InjOTCvsINF_NoTRMT'

#Create final significant comparisons table
final.sigtab.p <- sigtab.DNEG3.p.ji



######### Day -3 INF_NoTRMT vs INF_OralOTC ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "-3_INF_NoTRMT")
#20
sum(meta$Set == "-3_INF_OralOTC")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.p.De$Set
res.DNEG3.p.oi = results(FS9.DNEG3.p.De, contrast=c("Set", "-3_INF_OralOTC", "-3_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.DNEG3.p.De, contrast=c("Set","-3_INF_OralOTC", "-3_INF_NoTRMT")) 
sigtab.DNEG3.p.oi = res.DNEG3.p.oi[which(res.DNEG3.p.oi$padj < .05), ]
sigtab.DNEG3.p.oi = cbind(as(sigtab.DNEG3.p.oi, "data.frame"), as(tax_table(FS9.DNEG3.p)[rownames(sigtab.DNEG3.p.oi), ], "matrix"))
#Error in dimnames(x) <- dn : 
#length of 'dimnames' [1] not equal to array extent

sigtab.DNEG3.p.oi
#log2 fold change (MLE): Set -3_INF_OralOTC vs -3_INF_NoTRMT 
#Wald test p-value: Set -3_INF_OralOTC vs -3_INF_NoTRMT 
#DataFrame with 0 rows and 6 columns
#0 rows and 6 columns indicate there were no phyla that were found significantly differentially abundant between INF_OralOTC and INF_NoTRMT.
#I will continue to the next data set.

#SKIPPED THIS
format(sigtab.DNEG3.p.oi$padj, scientific = TRUE)
sigtab.DNEG3.p.oi$newp <- format(round(sigtab.DNEG3.p.oi$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.p.oi$Treatment <- ifelse(sigtab.DNEG3.p.oi$log2FoldChange >=0, "INF_OralOTC", "INF_NoTRMT")

#SKIPPED THIS
#Summarize sigtab.DNEG3.p.oi
sum.sigtab.DNEG3.p.oi <- summary(sigtab.DNEG3.p.oi)
sum.sigtab.DNEG3.p.oi

#SKIPPED THIS
#ggplot
deseq.DNEG3.p.oi <- ggplot(sigtab.DNEG3.p.oi, aes(x=reorder(rownames(sigtab.DNEG3.p.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p.oi), y=-1, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phyla in INF_OralOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_OralOTC", "INF_NoTRMT"), values = c('#999999', '#CC0066'))
  scale_fill_manual(values = c(INF_OralOTC='#999999', INF_NoTRMT='#CC0066'))
deseq.DNEG3.p.oi

#SKIPPED THIS
#Add OTU and comparisons columns
sigtab.DNEG3.p.oi
sigtab.DNEG3.p.oi$OTU <- rownames(sigtab.DNEG3.p.oi)
sigtab.DNEG3.p.oi
sigtab.DNEG3.p.oi$comp <- 'DNEG3_INF_OralOTCvsINF_NoTRMT'

#SKIPPED THIS
#Create final significant comparisons table
final.sigtab.p <- rbind(sigtab.DNEG3.p.oi, final.sigtab.p)



######### Day -3 NONINF_NoTRMT vs INF_NoTRMT  ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "-3_NONINF_NoTRMT")
#20
sum(meta$Set == "-3_INF_NoTRMT")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.p.De$Set
res.DNEG3.p.in = results(FS9.DNEG3.p.De, contrast=c("Set", "-3_INF_NoTRMT", "-3_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.DNEG3.p.De, contrast=c("Set","-3_INF_NoTRMT", "-3_NONINF_NoTRMT")) 
sigtab.DNEG3.p.in = res.DNEG3.p.in[which(res.DNEG3.p.in$padj < .05), ]
sigtab.DNEG3.p.in = cbind(as(sigtab.DNEG3.p.in, "data.frame"), as(tax_table(FS9.DNEG3.p)[rownames(sigtab.DNEG3.p.in), ], "matrix"))
format(sigtab.DNEG3.p.in$padj, scientific = TRUE)
sigtab.DNEG3.p.in$newp <- format(round(sigtab.DNEG3.p.in$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.p.in$Treatment <- ifelse(sigtab.DNEG3.p.in$log2FoldChange >=0, "INF_NoTRMT", "NONINF_NoTRMT")

#Summarize sigtab.DNEG3.p.in
sum.sigtab.DNEG3.p.in <- summary(sigtab.DNEG3.p.in)
sum.sigtab.DNEG3.p.in

#ggplot
deseq.DNEG3.p.in <- ggplot(sigtab.DNEG3.p.in, aes(x=reorder(rownames(sigtab.DNEG3.p.in), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p.in), y=1, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phyla in INF_NoTRMT Group Relative to NONINF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_NoTRMT", "NONINF_NoTRMT"), values = c('#CC0066', '#56B4E9'))
  scale_fill_manual(values = c(INF_NoTRMT='#CC0066', NONINF_NoTRMT='#56B4E9'))
deseq.DNEG3.p.in

#Add OTU and comparisons columns
sigtab.DNEG3.p.in
sigtab.DNEG3.p.in$OTU <- rownames(sigtab.DNEG3.p.in)
sigtab.DNEG3.p.in
sigtab.DNEG3.p.in$comp <- 'DNEG3_INF_NoTRMTvsNONINF_NoTRMT'

#Create final significant comparisons table
final.sigtab.p <- rbind(final.sigtab.p, sigtab.DNEG3.p.in)



######### Day -3 NONINF_NoTRMT vs INF_OralOTC  ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "-3_NONINF_NoTRMT")
#20
sum(meta$Set == "-3_INF_OralOTC")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.p.De$Set
res.DNEG3.p.on = results(FS9.DNEG3.p.De, contrast=c("Set", "-3_INF_OralOTC", "-3_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.DNEG3.p.De, contrast=c("Set","-3_INF_OralOTC", "-3_NONINF_NoTRMT")) 
sigtab.DNEG3.p.on = res.DNEG3.p.on[which(res.DNEG3.p.on$padj < .05), ]
sigtab.DNEG3.p.on = cbind(as(sigtab.DNEG3.p.on, "data.frame"), as(tax_table(FS9.DNEG3.p)[rownames(sigtab.DNEG3.p.on), ], "matrix"))
format(sigtab.DNEG3.p.on$padj, scientific = TRUE)
sigtab.DNEG3.p.on$newp <- format(round(sigtab.DNEG3.p.on$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.p.on$Treatment <- ifelse(sigtab.DNEG3.p.on$log2FoldChange >=0, "INF_OralOTC", "NONINF_NoTRMT")

#Summarize sigtab.DNEG3.p.on
sum.sigtab.DNEG3.p.on <- summary(sigtab.DNEG3.p.on)
sum.sigtab.DNEG3.p.on

#ggplot
deseq.DNEG3.p.on <- ggplot(sigtab.DNEG3.p.on, aes(x=reorder(rownames(sigtab.DNEG3.p.on), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p.on), y=1, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phyla in INF_OralOTC Group Relative to NONINF_NoTRMT in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_OralOTC", "NONINF_NoTRMT"), values = c('#999999', '#56B4E9'))
  scale_fill_manual(values = c(INF_OralOTC='#999999', NONINF_NoTRMT='#56B4E9'))
deseq.DNEG3.p.on

#Add OTU and comparisons columns
sigtab.DNEG3.p.on
sigtab.DNEG3.p.on$OTU <- rownames(sigtab.DNEG3.p.on)
sigtab.DNEG3.p.on
sigtab.DNEG3.p.on$comp <- 'DNEG3_INF_OralOTCvsNONINF_NoTRMT'

#Create final significant comparisons table
final.sigtab.p <- rbind(final.sigtab.p, sigtab.DNEG3.p.on)


##################################################### Day 7 ######################################################################

sample_data(FS9.phylum)

FS9.D7.p <- subset_samples(FS9.phylum, Day == '7')
sample_sums(FS9.D7.p)
colnames(otu_table(FS9.D7.p)) #check on all the sample names
FS9.D7.p <- prune_taxa(taxa_sums(FS9.D7.p) > 1, FS9.D7.p)
#if taxa_sums is >1, then it will print that out in FS9.D7.p object and not include anything with <1.
rowSums(FS9.D7.p@otu_table)

#Look at what Set is
sample_data(FS9.D7.p)

FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing

######### Day 7 INF_InjOTC vs INF_NoTRMT  ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC= J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "7_INF_InjOTC")
#26
sum(meta$Set == "7_INF_NoTRMT")
#20

#Extract results from a DESeq analysis, organize table
FS9.D7.p.De$Set
res.D7.p.ji = results(FS9.D7.p.De, contrast=c("Set", "7_INF_InjOTC", "7_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.D7.p.De, contrast=c("Set","7_INF_InjOTC", "7_INF_NoTRMT")) 
sigtab.D7.p.ji = res.D7.p.ji[which(res.D7.p.ji$padj < .05), ]
sigtab.D7.p.ji = cbind(as(sigtab.D7.p.ji, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p.ji), ], "matrix"))
format(sigtab.D7.p.ji$padj, scientific = TRUE)
sigtab.D7.p.ji$newp <- format(round(sigtab.D7.p.ji$padj, digits = 3), scientific = TRUE)
sigtab.D7.p.ji$Treatment <- ifelse(sigtab.D7.p.ji$log2FoldChange >=0, "INF_InjOTC", "INF_NoTRMT")

#Summarize sigtab.D7.p.ji
sum.sigtab.D7.p.ji <- summary(sigtab.D7.p.ji)
sum.sigtab.D7.p.ji

#ggplot
deseq.D7.p.ji <- ggplot(sigtab.D7.p.ji, aes(x=reorder(rownames(sigtab.D7.p.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.ji), y=-2, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phyla in INF_InjOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_InjOTC", "INF_NoTRMT"), values = c('#E69F00', '#CC0066'))
  scale_fill_manual(values = c(INF_InjOTC='#E69F00', INF_NoTRMT='#CC0066'))
deseq.D7.p.ji

#Add OTU and comparisons columns
sigtab.D7.p.ji
sigtab.D7.p.ji$OTU <- rownames(sigtab.D7.p.ji)
sigtab.D7.p.ji
sigtab.D7.p.ji$comp <- 'D7_INF_InjOTCvsINF_NoTRMT'

#Create final significant comparisons table
final.sigtab.p <- rbind(sigtab.D7.p.ji, final.sigtab.p)



######### Day 7 INF_OralOTC vs INF_NoTRMT ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "7_INF_NoTRMT")
#20
sum(meta$Set == "7_INF_OralOTC")
#17

#Extract results from a DESeq analysis, organize table
FS9.D7.p.De$Set
res.D7.p.oi = results(FS9.D7.p.De, contrast=c("Set", "7_INF_OralOTC", "7_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.D7.p.De, contrast=c("Set","7_INF_OralOTC", "7_INF_NoTRMT")) 
sigtab.D7.p.oi = res.D7.p.oi[which(res.D7.p.oi$padj < .05), ]
sigtab.D7.p.oi = cbind(as(sigtab.D7.p.oi, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p.oi), ], "matrix"))
format(sigtab.D7.p.oi$padj, scientific = TRUE)
sigtab.D7.p.oi$newp <- format(round(sigtab.D7.p.oi$padj, digits = 3), scientific = TRUE)
sigtab.D7.p.oi$Treatment <- ifelse(sigtab.D7.p.oi$log2FoldChange >=0, "INF_OralOTC", "INF_NoTRMT")

#Summarize sigtab.D7.p.oi
sum.sigtab.D7.p.oi <- summary(sigtab.D7.p.oi)
sum.sigtab.D7.p.oi

#ggplot
deseq.D7.p.oi <- ggplot(sigtab.D7.p.oi, aes(x=reorder(rownames(sigtab.D7.p.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.oi), y=-3, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phyla in INF_OralOTC Group Relative to INF_NoTRMT in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_OralOTC", "INF_NoTRMT"), values = c('#999999', '#CC0066'))
  scale_fill_manual(values = c(INF_OralOTC='#999999', INF_NoTRMT='#CC0066'))
deseq.D7.p.oi

#Add OTU and comparisons columns
sigtab.D7.p.oi
sigtab.D7.p.oi$OTU <- rownames(sigtab.D7.p.oi)
sigtab.D7.p.oi
sigtab.D7.p.oi$comp <- 'D7_INF_OralOTCvsINF_NoTRMT'

#Create final significant comparisons table
final.sigtab.p <- rbind(final.sigtab.p, sigtab.D7.p.oi)




######### Day 7 INF_OralOTC vs INF_InjOTC ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC = J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "7_INF_InjOTC")
#26
sum(meta$Set == "7_INF_OralOTC")
#17

#Extract results from a DESeq analysis, organize table
FS9.D7.p.De$Set
res.D7.p.oj = results(FS9.D7.p.De, contrast=c("Set", "7_INF_OralOTC", "7_INF_InjOTC"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.D7.p.De, contrast=c("Set","7_INF_OralOTC", "7_INF_InjOTC")) 
sigtab.D7.p.oj = res.D7.p.oj[which(res.D7.p.oj$padj < .05), ]
sigtab.D7.p.oj = cbind(as(sigtab.D7.p.oj, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p.oj), ], "matrix"))
format(sigtab.D7.p.oj$padj, scientific = TRUE)
sigtab.D7.p.oj$newp <- format(round(sigtab.D7.p.oj$padj, digits = 3), scientific = TRUE)
sigtab.D7.p.oj$Treatment <- ifelse(sigtab.D7.p.oj$log2FoldChange >=0, "INF_OralOTC", "INF_InjOTC")

#Summarize sigtab.D7.p.oj
sum.sigtab.D7.p.oj <- summary(sigtab.D7.p.oj)
sum.sigtab.D7.p.oj

#ggplot
deseq.D7.p.oj <- ggplot(sigtab.D7.p.oj, aes(x=reorder(rownames(sigtab.D7.p.oj), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.oj), y=1, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phyla in INF_OralOTC Group Relative to INF_InjOTC in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_OralOTC", "INF_InjOTC"), values = c('#999999', '#E69F00'))
  scale_fill_manual(values = c(INF_OralOTC='#999999', INF_InjOTC='#E69F00'))
deseq.D7.p.oj

#Add OTU and comparisons columns
sigtab.D7.p.oj
sigtab.D7.p.oj$OTU <- rownames(sigtab.D7.p.oj)
sigtab.D7.p.oj
sigtab.D7.p.oj$comp <- 'D7_INF_OralOTCvsINF_InjOTC'

#Create final significant comparisons table
final.sigtab.p <- rbind(final.sigtab.p, sigtab.D7.p.oj)



######### Day 7 INF_InjOTC vs NONINF_NoTRMT  ###################

#NONINF_NoTRMT = N
#INF_NoTRMT = I
#INF_InjOTC= J
#INF_OralOTC = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "7_INF_InjOTC")
#26
sum(meta$Set == "7_NONINF_NoTRMT")
#17

#Extract results from a DESeq analysis, organize table
FS9.D7.p.De$Set
res.D7.p.jn = results(FS9.D7.p.De, contrast=c("Set", "7_INF_InjOTC", "7_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.D7.p.De, contrast=c("Set","7_INF_InjOTC", "7_NONINF_NoTRMT")) 
sigtab.D7.p.jn = res.D7.p.jn[which(res.D7.p.jn$padj < .05), ]
sigtab.D7.p.jn = cbind(as(sigtab.D7.p.jn, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p.jn), ], "matrix"))
format(sigtab.D7.p.jn$padj, scientific = TRUE)
sigtab.D7.p.jn$newp <- format(round(sigtab.D7.p.jn$padj, digits = 3), scientific = TRUE)
sigtab.D7.p.jn$Treatment <- ifelse(sigtab.D7.p.jn$log2FoldChange >=0, "INF_InjOTC", "NONINF_NoTRMT")

#Summarize sigtab.D7.p.jn
sum.sigtab.D7.p.jn <- summary(sigtab.D7.p.jn)
sum.sigtab.D7.p.jn

#ggplot
deseq.D7.p.jn <- ggplot(sigtab.D7.p.jn, aes(x=reorder(rownames(sigtab.D7.p.jn), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.jn), y=-2, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phyla in INF_InjOTC Group Relative to NONINF_NoTRMT in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  #scale_fill_manual(labels = c("INF_InjOTC", "NONINF_NoTRMT"), values = c('#E69F00', '#56B4E9'))
  scale_fill_manual(values = c(INF_InjOTC='#E69F00', NONINF_NoTRMT='#56B4E9'))
deseq.D7.p.jn

#Add OTU and comparisons columns
sigtab.D7.p.jn
sigtab.D7.p.jn$OTU <- rownames(sigtab.D7.p.jn)
sigtab.D7.p.jn
sigtab.D7.p.jn$comp <- 'D7_INF_InjOTCvsNONINF_NoTRMT'

#Create final significant comparisons table
final.sigtab.p <- rbind(sigtab.D7.p.jn, final.sigtab.p)



#######################################################################################################

#write csv
write.csv(final.sigtab.p, file= "FS9_FinalDiffAbund_SignificantDays_Phylum.csv")

#######################################################################################################





######### Plots of Differentially Abundant Nasal Genera Combined for Each Pairwise Comparison of IAV and Control Groups########
library("ggsci")

#IAV and Control Log2fold Plot Part A
final_ic <- sigtab.DNEG3.ji
final_ic <- rbind(final_ic, sigtab.DNEG12.ic, sigtab.DNEG6.ic, sigtab.D1.ic, sigtab.D3.ic, sigtab.D7.ic, sigtab.D10.ic, 
                  sigtab.D14.ic, sigtab.D21.ic, sigtab.D28.ic, sigtab.D36.ic, sigtab.D42.ic)
final_ic$Family_Genus <- paste(final_ic$Family, final_ic$Genus) #create new column with Family_Genus
final_ic$comp
class(final_ic)
final_ic$comp <- factor(final_ic$comp, levels=c("DNEG12_nasal_IAVvsControl", "DNEG6_nasal_IAVvsControl", "D0_nasal_IAVvsControl",
                                                "D1_nasal_IAVvsControl", "D3_nasal_IAVvsControl", "D7_nasal_IAVvsControl",
                                                "D10_nasal_IAVvsControl", "D14_nasal_IAVvsControl", "D21_nasal_IAVvsControl",
                                                "D28_nasal_IAVvsControl", "D36_nasal_IAVvsControl", "D42_nasal_IAVvsControl"))
levels(final_ic$comp)
ic_plot <- ggplot(final_ic, aes(x=Family_Genus, log2FoldChange, fill = comp)) +
  geom_bar(stat='identity') +
  labs(x="Family Genus", y = "Total log2 Fold Change") +
  theme(axis.text.x=element_text(color = 'black', size = 18),
        axis.text.y=element_text(color = 'black', size=15), 
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20))+ 
  coord_flip() +
  scale_fill_igv(name = "comp") +
  ggtitle('Differentially Abundant Nasal Families and Genera between Influenza A Virus and Control Groups') + 
  theme(plot.title = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=20))
ic_plot <- ic_plot + guides(fill=guide_legend(title="Day and treatment group"))
ic_plot

write.csv(final_ic, file= "SRD129_FinalDiffAbundNasalGenus_IAVcontrol.csv")

#Modified SRD129_FinalDiffAbundNasalGenus_IAVcontrol.csv in a spreadsheet editor by removing all genera except for Actinobacillus, Moraxella, Neisseriaceae_unclassified, Prevotellaceae_NK3B31_group,
#Prevotellaceae_unclassified, Staphylococcus, Streptococcus.
#Saved modified file as SRD129_FinalDiffAbundNasalGenus_IAVcontrolSelectList.csv

#IAV and control Log2fold Plot Part B
ic <- read.csv('SRD129_FinalDiffAbundNasalGenus_IAVcontrolSelectList.csv', header = TRUE, sep = ",")
head(ic[,1:10])
colnames(ic)
ic$DayComp <- sub('_[A-Za-z]+', '\\2', ic$comp)
unique(ic$DayComp)
ic$Day <- sub('_[A-Za-z]+', '\\2', ic$DayComp)
unique(ic$Day) #"D3"  "D7"  "D10" "D14" "D21" "D28"
ic$Day = factor(ic$Day, levels=c("D3","D7", "D10","D14", "D21", "D28"))
levels(ic$Day) #"D3"  "D7"  "D10" "D14" "D21" "D28"
(ic_logfoldplot <- ggplot(data=ic, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 3, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
          axis.text.y = element_text(size=13), 
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15),
          legend.text=element_text(size=15), 
          legend.title=element_text(size=15)) +
    scale_fill_manual(labels = c("Control", "IAV"), values = c('#F8766D', '#00BFC4')))
ic_logfoldplot <- ic_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
ic_logfoldplot

#Save 'ic_logfoldplot' as a .tiff for publication, 500dpi
ggsave("Figure_5.tiff", plot=ic_logfoldplot, width = 15, height = 7, dpi = 500, units =c("in"))

#Figure for Sam (All days except negative days, IAV Control complete list of diff. abund genera)
id <- read.csv('SRD129_FinalDiffAbundNasalGenus_IAVcontrol_NoDNEG.csv', header = TRUE, sep = ",")
head(id[,1:10])
colnames(id)
id$DayComp <- sub('_[A-Za-z]+', '\\2', id$comp)
unique(id$DayComp)
id$Day <- sub('_[A-Za-z]+', '\\2', id$DayComp)
unique(id$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
id$Day = factor(id$Day, levels=c("D0", "D1", "D3", "D7", "D10","D14", "D21", "D28", "D36", "D42"))
levels(id$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
(id_logfoldplot <- ggplot(data=id, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 10, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
          axis.text.y = element_text(size=13), 
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15),
          legend.text=element_text(size=15), 
          legend.title=element_text(size=15)) +
    scale_fill_manual(labels = c("Control", "IAV"), values = c('#F8766D', '#00BFC4')))
id_logfoldplot <- id_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
id_logfoldplot
unique(id$Genus)

#Figure for Sam, narrow down complete list to 46 genera (All days except negative days)
ie <- read.csv('SRD129_FinalDiffAbundNasalGenus_IAVcontrol_NoDNEG_46Genera.csv', header = TRUE, sep = ",")
head(ie[,1:10])
colnames(ie)
ie$DayComp <- sub('_[A-Za-z]+', '\\2', ie$comp)
unique(ie$DayComp)
ie$Day <- sub('_[A-Za-z]+', '\\2', ie$DayComp)
unique(ie$Day) # "D7"  "D14" "D3"  "D10" "D21" "D36" "D42" "D0"  "D28" "D1" 
ie$Day = factor(ie$Day, levels=c("D0", "D1", "D3", "D7", "D10","D14", "D21", "D28", "D36", "D42"))
levels(ie$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
(ie_logfoldplot <- ggplot(data=ie, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 8, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
          axis.text.y = element_text(size=13), 
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15),
          legend.text=element_text(size=15), 
          legend.title=element_text(size=15)) +
    scale_fill_manual(labels = c("Control", "IAV"), values = c('#F8766D', '#00BFC4')))
ie_logfoldplot <- ie_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
ie_logfoldplot
unique(ie$Genus)

#Figure for Sam, narrow down the 46 genera to only those that aren't diff abund on day 0 (All days except negative days)
iF <- read.csv('SRD129_FinalDiffAbundNasalGenus_IAVcontrol_NoDNEG_46Genera-Day0Genera.csv', header = TRUE, sep = ",")
head(iF[,1:10])
colnames(iF)
iF$DayComp <- sub('_[A-Za-z]+', '\\2', iF$comp)
unique(iF$DayComp)
iF$Day <- sub('_[A-Za-z]+', '\\2', iF$DayComp)
unique(iF$Day) # "D7"  "D14" "D3"  "D10" "D21" "D36" "D42" "D0"  "D28" "D1" 
iF$Day = factor(iF$Day, levels=c("D0", "D1", "D3", "D7", "D10","D14", "D21", "D28", "D36", "D42"))
levels(iF$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
(iF_logfoldplot <- ggplot(data=iF, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 8, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
          axis.text.y = element_text(size=13), 
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15),
          legend.text=element_text(size=15), 
          legend.title=element_text(size=15)) +
    scale_fill_manual(labels = c("Control", "IAV"), values = c('#F8766D', '#00BFC4')))
iF_logfoldplot <- iF_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
iF_logfoldplot
unique(iF$Genus)
