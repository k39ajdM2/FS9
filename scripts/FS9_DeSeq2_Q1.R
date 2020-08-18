#####################################################################################################
#FS9 DESeq2 - Noninfected vs Infected
#Kathy Mou

#Purpose: This code uses DESeq2 package to identify fecal microbial genera that were differentially 
#abundant between the two groups on the two days sampled

#Load library packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("phyloseq")
library(phyloseq)
library(ggplot2)
library("wesanderson")
library(plotly)
library(gapminder)
library("ggsci")
BiocManager::install("apeglm")
library("apeglm")

sessionInfo()
#R version 3.6.0 (2019-04-26)

########################################################################################################

#Import files
otu <- import_mothur(mothur_shared_file = './data/stability.outdoubletons.abund.opti_mcc.shared') #use unrarified data
taxo <- import_mothur(mothur_constaxonomy_file = './data/stability.outdoubletons.abund.opti_mcc.0.03.cons.taxonomy')
meta <- read.table(file = './data/FS9_metadata_NONINFvINF.csv', sep = ',', header = TRUE)

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
#otu_table()   OTU Table:         [ 2039 taxa and 112 samples ]
#sample_data() Sample Data:       [ 112 samples by 5 sample variables ]
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
FS9.phylum <- tax_glom(FS9, taxrank = "Phylum")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 

######################################## PRIMARY COMPARISONS TO MAKE ############################################################

# NMDS plot showed that disperion is different between days, so I subsetted by day
# Though there were no significant differences between groups on all days sampled, I will still look at
# taxons to see if there were consistent differences in abundance between groups across days

# Comparisons to make:
# Day -3 
# INF vs NONINF

# Day 0 
# INF vs NONINF

######################################################### Day -3 #########################################################

unique(sample_data(FS9.order)$Set)

FS9.DNEG3 <- subset_samples(FS9.order, Day == 'DNEG3')
sample_sums(FS9.DNEG3)
colnames(otu_table(FS9.DNEG3)) #check on all the sample names
FS9.DNEG3 <- prune_taxa(taxa_sums(FS9.DNEG3) > 1, FS9.DNEG3)
#if taxa_sums is >1, then it will print that out in FS9.DNEG3 object and not include anything with <1.
rowSums(FS9.DNEG3@otu_table)

#Look at what Set is
unique(sample_data(FS9.DNEG3)$Set)

FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")

######### Day -3 INF vs NONINF ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INF")
#60
sum(meta$Set == "DNEG3_NONINF")
#20

resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                            "Set_DNEG3_NONINF_vs_DNEG3_INF"

#re-level your factor and re-run DESeq2
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c('DNEG3_NONINF',"DNEG3_INF"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                "Set_DNEG3_INF_vs_DNEG3_NONINF"
res.DNEG3 = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INF_vs_DNEG3_NONINF", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.DNEG3 = res.DNEG3[which(res.DNEG3$padj < .05), ]
sigtab.DNEG3 = cbind(as(sigtab.DNEG3, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3), ], "matrix"))
sigtab.DNEG3$newp <- format(round(sigtab.DNEG3$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3$Treatment <- ifelse(sigtab.DNEG3$log2FoldChange >=0, "INF", "NONINF")

deseq.DNEG3 <- 
  ggplot(sigtab.DNEG3, aes(x=reorder(rownames(sigtab.DNEG3), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3), y=1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INF Group Relative to NONINF\n in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INF='#CC0066', NONINF='#56B4E9'))
deseq.DNEG3

#Add OTU and comparisons columns
sigtab.DNEG3
sigtab.DNEG3$OTU <- rownames(sigtab.DNEG3)
sigtab.DNEG3
sigtab.DNEG3$comp <- 'DNEG3_INF_vs_NONINF'

#Create final significant comparisons table
final.sigtab <- sigtab.DNEG3

##################################################### Day 0 ######################################################################

FS9.D0 <- subset_samples(FS9.order, Day == 'D0')
sample_sums(FS9.D0)
colnames(otu_table(FS9.D0)) #check on all the sample names
FS9.D0 <- prune_taxa(taxa_sums(FS9.D0) > 1, FS9.D0)
#if taxa_sums is >1, then it will print that out in FS9.D0 object and not include anything with <1.
rowSums(FS9.D0@otu_table)

#Look at what Set is
unique(sample_data(FS9.D0)$Set)

FS9.D0.De <- phyloseq_to_deseq2(FS9.D0, ~ Set)
FS9.D0.De <- DESeq(FS9.D0.De, test = "Wald", fitType = "parametric")

######### Day 0 INF vs NONINF ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_INF")
#26
sum(meta$Set == "D0_NONINF")
#9

resultsNames(FS9.D0.De)
#[1] [1] "Intercept"            "Set_D0_NONINF_vs_D0_INF"
sample_data(FS9.D0)$Set <- factor(sample_data(FS9.D0)$Set,
                                     levels =c('D0_NONINF',"D0_INF"))
FS9.D0.De <- phyloseq_to_deseq2(FS9.D0, ~ Set)
FS9.D0.De <- DESeq(FS9.D0.De, test = "Wald", fitType = "parametric")
FS9.D0.De$Set
resultsNames(FS9.D0.De)
#[1] "Intercept"            "Set_D0_INF_vs_D0_NONINF"
res.D0 = lfcShrink(FS9.D0.De, coef = "Set_D0_INF_vs_D0_NONINF", type = 'apeglm')
sigtab.D0 = res.D0[which(res.D0$padj < .05), ]
sigtab.D0 = cbind(as(sigtab.D0, "data.frame"), as(tax_table(FS9.D0)[rownames(sigtab.D0), ], "matrix"))
sigtab.D0$newp <- format(round(sigtab.D0$padj, digits = 3), scientific = TRUE)
sigtab.D0$Treatment <- ifelse(sigtab.D0$log2FoldChange >=0, "INF", "NONINF")

deseq.D0 <- 
  ggplot(sigtab.D0, aes(x=reorder(rownames(sigtab.D0), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0), y=0, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INF Group Relative to NONINF\n in Fecal Microbiota on Day 0')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INF='#CC0066', NONINF='#56B4E9'))
deseq.D0

#Add OTU and comparisons columns
sigtab.D0
sigtab.D0$OTU <- rownames(sigtab.D0)
sigtab.D0
sigtab.D0$comp <- 'D0_INF_vs_NONINF'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.D0, final.sigtab)

#Filtering out OTUs that have log2-fold values less than 0.25
final.sigtab2 <- final.sigtab %>% 
  filter(log2FoldChange > 0.25) 
deseq.D0.DNEG3 <- 
  ggplot(final.sigtab2, aes(x=reorder(rownames(final.sigtab2), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(final.sigtab2), y=1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INF Group Relative to NONINF\n in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INF='#CC0066', NONINF='#56B4E9'))
deseq.D0.DNEG3

#######################################################################################################

#write csv that includes all differentially abundant OTUs regardless of log2-fold values
write.csv(final.sigtab, file= "FS9_FinalDiffAbund_Order_OutDoubletons_Q1.csv")

#######################################################################################################






##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (PHYLUM LEVEL) BETWEEN TREATMENTS ###################################

######################################################### Day -3 #########################################################

sample_data(FS9.phylum)

FS9.DNEG3.p <- subset_samples(FS9.phylum, Day == 'DNEG3')
sample_sums(FS9.DNEG3.p)
colnames(otu_table(FS9.DNEG3.p)) #check on all the sample names
FS9.DNEG3.p <- prune_taxa(taxa_sums(FS9.DNEG3.p) > 1, FS9.DNEG3.p)
#if taxa_sums is >1, then it will print that out in FS9.DNEG3.p object and not include anything with <1.
rowSums(FS9.DNEG3.p@otu_table)

#Look at what Set is
sample_data(FS9.DNEG3.p)

######### 1. Day -3 INFinject vs INFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INF")
#60
sum(meta$Set == "DNEG3_NONINF")
#20

#re-level your factor and re-run DESeq2
sample_data(FS9.DNEG3.p)$Set <- factor(sample_data(FS9.DNEG3.p)$Set,
                                     levels =c('DNEG3_NONINF',"DNEG3_INF"))
FS9.DNEG3.p.De <- phyloseq_to_deseq2(FS9.DNEG3.p, ~ Set)
FS9.DNEG3.p.De <- DESeq(FS9.DNEG3.p.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.p.De$Set
resultsNames(FS9.DNEG3.p.De)
#[1] "Intercept"                   "Set_DNEG3_INF_vs_DNEG3_NONINF"  
res.DNEG3.p = lfcShrink(FS9.DNEG3.p.De, coef = "Set_DNEG3_INF_vs_DNEG3_NONINF", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (DNEG3_INF)
# negative log2foldchanges are associated with the second group from this line (DNEG3_NONINF)
sigtab.DNEG3.p = res.DNEG3.p[which(res.DNEG3.p$padj < .05), ]
sigtab.DNEG3.p = cbind(as(sigtab.DNEG3.p, "data.frame"), as(tax_table(FS9.DNEG3.p)[rownames(sigtab.DNEG3.p), ], "matrix"))
sigtab.DNEG3.p$newp <- format(round(sigtab.DNEG3.p$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.p$Treatment <- ifelse(sigtab.DNEG3.p$log2FoldChange >=0, "INF", "NONINF")
head(sigtab.DNEG3.p)

deseq.DNEG3.p <- 
  ggplot(sigtab.DNEG3.p, aes(x=reorder(rownames(sigtab.DNEG3.p), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p), y=0, label = paste(Phylum)), size=5, fontface = 'italic')+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phylum\n in INF Group Relative to NONINF\n in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INF='#CC0066', NONINF='#56B4E9'))
deseq.DNEG3.p

#Add OTU and comparisons columns
sigtab.DNEG3.p
sigtab.DNEG3.p$OTU <- rownames(sigtab.DNEG3.p)
sigtab.DNEG3.p
sigtab.DNEG3.p$comp <- 'DNEG3_INF_vs_NONINF'

#Create final significant comparisons table
final.sigtab.phylum <- sigtab.DNEG3.p



##################################################### Day 0 ######################################################################

FS9.D0.p <- subset_samples(FS9.phylum, Day == 'D0')
sample_sums(FS9.D0.p)
colnames(otu_table(FS9.D0.p)) #check on all the sample names
FS9.D0.p <- prune_taxa(taxa_sums(FS9.D0.p) > 1, FS9.D0.p)
#if taxa_sums is >1, then it will print that out in FS9.D0 object and not include anything with <1.
rowSums(FS9.D0.p@otu_table)

#Look at what Set is
unique(sample_data(FS9.D0.p)$Set)

######### 1. Day 0 INFinject vs INFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_INF")
#26
sum(meta$Set == "D0_NONINF")
#9

sample_data(FS9.D0.p)$Set <- factor(sample_data(FS9.D0.p)$Set,
                                  levels =c('D0_NONINF',"D0_INF"))
FS9.D0.p.De <- phyloseq_to_deseq2(FS9.D0.p, ~ Set)
FS9.D0.p.De <- DESeq(FS9.D0.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D0.p.De)
#[1] "Intercept"     "Set_D0_INF_vs_D0_NONINF"
res.D0.p = lfcShrink(FS9.D0.p.De, coef = "Set_D0_INF_vs_D0_NONINF", type = 'apeglm')
sigtab.D0.p = res.D0.p[which(res.D0.p$padj < .05), ]
sigtab.D0.p = cbind(as(sigtab.D0.p, "data.frame"), as(tax_table(FS9.D0.p)[rownames(sigtab.D0.p), ], "matrix"))
sigtab.D0.p$newp <- format(round(sigtab.D0.p$padj, digits = 3), scientific = TRUE)
sigtab.D0.p$Treatment <- ifelse(sigtab.D0.p$log2FoldChange >=0, "INF", "NONINF")
head(sigtab.D0.p) 

#CONTINUE HERE!



#######################################################################################################

#write csv
write.csv(final.sigtab.phylum, file= "FS9_FinalDiffAbund_Phylum_OutDoubletons_Q1.csv")

#######################################################################################################