#####################################################################################################
#FS9 DESeq2 - Noninfected (NONINFnm) vs Infected (INFnm) day -3; INFnm vs INFinject vs INFfeed
#Kathy Mou

#Purpose: This code uses DESeq2 package to identify fecal microbial genera that were differentially 
#abundant between the two groups NONINFnm and INFnm on all days and only days 4 and 7 at the order level

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
#R version 4.0.2 (2020-06-22)

########################################################################################################

#Import files
otu <- import_mothur(mothur_shared_file = './data/stability.outdoubletons.abund.opti_mcc.shared') #use unrarified data
taxo <- import_mothur(mothur_constaxonomy_file = './data/stability.outdoubletons.abund.opti_mcc.0.03.cons.taxonomy')
meta <- read.table(file = './data/FS9_metadata_NONINFvINF_alldays.csv', sep = ',', header = TRUE)

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
#otu_table()   OTU Table:         [ 2039 taxa and 85 samples ]
#sample_data() Sample Data:       [ 85 samples by 5 sample variables ]
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
FS9.genus <- tax_glom(FS9, taxrank = "Genus")

##################################################################################################################################

##################################### Q1: SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (ORDER LEVEL) BETWEEN TREATMENTS ###################################

##################################### Day -3, 0, 4, 7 ##########################################################################


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

######### Day -3 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFnm")
#20
sum(meta$Set == "DNEG3_NONINFnm")
#20

resultsNames(FS9.DNEG3.De)
#[1] "Intercept"               "Set_DNEG3_NONINFnm_vs_DNEG3_INFnm"

#re-level your factor and re-run DESeq2
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c('DNEG3_NONINFnm',"DNEG3_INFnm"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"
res.DNEG3 = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.DNEG3 = res.DNEG3[which(res.DNEG3$padj < .05), ]
sigtab.DNEG3 = cbind(as(sigtab.DNEG3, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3), ], "matrix"))
sigtab.DNEG3$newp <- format(round(sigtab.DNEG3$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3$Treatment <- ifelse(sigtab.DNEG3$log2FoldChange >=0, "INFnm", "NONINFnm")

deseq.DNEG3 <- 
  ggplot(sigtab.DNEG3, aes(x=reorder(rownames(sigtab.DNEG3), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3), y=1, label = paste(Phylum, Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.DNEG3

#Add OTU and comparisons columns
sigtab.DNEG3
sigtab.DNEG3$OTU <- rownames(sigtab.DNEG3)
sigtab.DNEG3
sigtab.DNEG3$comp <- 'DNEG3_INFnm_vs_NONINFnm'

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

######### Day 0 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_INFnm")
#8
sum(meta$Set == "D0_NONINFnm")
#9

resultsNames(FS9.D0.De)
#[1] [1] "Intercept"           "Set_D0_NONINFnm_vs_D0_INFnm"
sample_data(FS9.D0)$Set <- factor(sample_data(FS9.D0)$Set,
                                  levels =c('D0_NONINFnm',"D0_INFnm"))
FS9.D0.De <- phyloseq_to_deseq2(FS9.D0, ~ Set)
FS9.D0.De <- DESeq(FS9.D0.De, test = "Wald", fitType = "parametric")
FS9.D0.De$Set
resultsNames(FS9.D0.De)
#[1] "Intercept"             "Set_D0_INFnm_vs_D0_NONINFnm"
res.D0 = lfcShrink(FS9.D0.De, coef = "Set_D0_INFnm_vs_D0_NONINFnm", type = 'apeglm')
sigtab.D0 = res.D0[which(res.D0$padj < .05), ]
sigtab.D0 = cbind(as(sigtab.D0, "data.frame"), as(tax_table(FS9.D0)[rownames(sigtab.D0), ], "matrix"))
sigtab.D0$newp <- format(round(sigtab.D0$padj, digits = 3), scientific = TRUE)
sigtab.D0$Treatment <- ifelse(sigtab.D0$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D0)
#DataFrame with 0 rows and 7 columns. Will continue to next data set



######################################################### Day 4 #########################################################

unique(sample_data(FS9.order)$Set)

FS9.D4 <- subset_samples(FS9.order, Day == 'D4')
sample_sums(FS9.D4)
colnames(otu_table(FS9.D4)) #check on all the sample names
FS9.D4 <- prune_taxa(taxa_sums(FS9.D4) > 1, FS9.D4)
#if taxa_sums is >1, then it will print that out in FS9.D4 object and not include anything with <1.
rowSums(FS9.D4@otu_table)

#Look at what Set is
unique(sample_data(FS9.D4)$Set)

FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")

######### Day 4 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_INFnm")
#6
sum(meta$Set == "D4_NONINFnm")
#5

resultsNames(FS9.D4.De)
#[1] "Intercept"             "Set_D4_NONINFnm_vs_D4_INFnm"

#re-level your factor and re-run DESeq2
sample_data(FS9.D4)$Set <- factor(sample_data(FS9.D4)$Set,
                                  levels =c('D4_NONINFnm',"D4_INFnm"))
FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")
FS9.D4.De$Set
resultsNames(FS9.D4.De)
#[1] "Intercept"              "Set_D4_INFnm_vs_D4_NONINFnm"
res.D4 = lfcShrink(FS9.D4.De, coef = "Set_D4_INFnm_vs_D4_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.D4 = res.D4[which(res.D4$padj < .05), ]
sigtab.D4 = cbind(as(sigtab.D4, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4), ], "matrix"))
sigtab.D4$newp <- format(round(sigtab.D4$padj, digits = 3), scientific = TRUE)
sigtab.D4$Treatment <- ifelse(sigtab.D4$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D4)
#DataFrame with 0 rows and 7 columns. Will continue to next data set


######################################################### Day 7 #########################################################

unique(sample_data(FS9.order)$Set)

FS9.D7 <- subset_samples(FS9.order, Day == 'D7')
sample_sums(FS9.D7)
colnames(otu_table(FS9.D7)) #check on all the sample names
FS9.D7 <- prune_taxa(taxa_sums(FS9.D7) > 1, FS9.D7)
#if taxa_sums is >1, then it will print that out in FS9.D7 object and not include anything with <1.
rowSums(FS9.D7@otu_table)

#Look at what Set is
unique(sample_data(FS9.D7)$Set)

FS9.D7.De <- phyloseq_to_deseq2(FS9.D7, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.D7.De <- DESeq(FS9.D7.De, test = "Wald", fitType = "parametric")

######### Day 7 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFnm")
#14
sum(meta$Set == "D7_NONINFnm")
#12

sample_data(FS9.D7)$Set <- factor(sample_data(FS9.D7)$Set,
                                  levels =c('D7_NONINFnm',"D7_INFnm"))
FS9.D7.De <- phyloseq_to_deseq2(FS9.D7, ~ Set)
FS9.D7.De <- DESeq(FS9.D7.De, test = "Wald", fitType = "parametric")
FS9.D7.De$Set
resultsNames(FS9.D7.De)
#[1] "Intercept"              "Set_D7_INFnm_vs_D7_NONINFnm"
res.D7 = lfcShrink(FS9.D7.De, coef = "Set_D7_INFnm_vs_D7_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.D7 = res.D7[which(res.D7$padj < .05), ]
sigtab.D7 = cbind(as(sigtab.D7, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7), ], "matrix"))
sigtab.D7$newp <- format(round(sigtab.D7$padj, digits = 3), scientific = TRUE)
sigtab.D7$Treatment <- ifelse(sigtab.D7$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D7)

deseq.D7 <- 
  ggplot(sigtab.D7, aes(x=reorder(rownames(sigtab.D7), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7), y=1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.D7

#Add OTU and comparisons columns
sigtab.D7
sigtab.D7$OTU <- rownames(sigtab.D7)
sigtab.D7
sigtab.D7$comp <- 'D7_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D7)


#######################################################################################################

#write csv that includes all differentially abundant OTUs regardless of log2-fold values
write.csv(final.sigtab, file= "FS9_FinalDiffAbund_Order_OutDoubletons_Q1_AllDays.csv")

#######################################################################################################




##################################### Q1: SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (GENUS LEVEL) BETWEEN TREATMENTS ###################################

##################################### Day -3, 0, 4, 7 ##########################################################################


######################################################### Day -3 #########################################################

unique(sample_data(FS9.genus)$Set)

FS9.DNEG3.g <- subset_samples(FS9.genus, Day == 'DNEG3')
sample_sums(FS9.DNEG3.g)
colnames(otu_table(FS9.DNEG3.g)) #check on all the sample names
FS9.DNEG3.g <- prune_taxa(taxa_sums(FS9.DNEG3.g) > 1, FS9.DNEG3.g)
#if taxa_sums is >1, then it will print that out in FS9.DNEG3.g object and not include anything with <1.
rowSums(FS9.DNEG3.g@otu_table)

#Look at what Set is
unique(sample_data(FS9.DNEG3.g)$Set)

FS9.DNEG3.g.De <- phyloseq_to_deseq2(FS9.DNEG3.g, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.DNEG3.g.De <- DESeq(FS9.DNEG3.g.De, test = "Wald", fitType = "parametric")

######### Day -3 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFnm")
#20
sum(meta$Set == "DNEG3_NONINFnm")
#20

resultsNames(FS9.DNEG3.g.De)
#[1] "Intercept"               "Set_DNEG3_NONINFnm_vs_DNEG3_INFnm"

#re-level your factor and re-run DESeq2
sample_data(FS9.DNEG3.g)$Set <- factor(sample_data(FS9.DNEG3.g)$Set,
                                     levels =c('DNEG3_NONINFnm',"DNEG3_INFnm"))
FS9.DNEG3.g.De <- phyloseq_to_deseq2(FS9.DNEG3.g, ~ Set)
FS9.DNEG3.g.De <- DESeq(FS9.DNEG3.g.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.g.De$Set
resultsNames(FS9.DNEG3.g.De)
#[1] "Intercept"                "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"
res.DNEG3.g = lfcShrink(FS9.DNEG3.g.De, coef = "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.DNEG3.g = res.DNEG3.g[which(res.DNEG3.g$padj < .05), ]
sigtab.DNEG3.g = cbind(as(sigtab.DNEG3.g, "data.frame"), as(tax_table(FS9.DNEG3.g)[rownames(sigtab.DNEG3.g), ], "matrix"))
sigtab.DNEG3.g$newp <- format(round(sigtab.DNEG3.g$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.g$Treatment <- ifelse(sigtab.DNEG3.g$log2FoldChange >=0, "INFnm", "NONINFnm")

deseq.DNEG3.g <- 
  ggplot(sigtab.DNEG3.g, aes(x=reorder(rownames(sigtab.DNEG3.g), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.g), y=1, label = paste(Order,Genus, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Order Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.DNEG3.g

#Add OTU and comparisons columns
sigtab.DNEG3.g
sigtab.DNEG3.g$OTU <- rownames(sigtab.DNEG3.g)
sigtab.DNEG3.g
sigtab.DNEG3.g$comp <- 'DNEG3_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab.g <- sigtab.DNEG3.g

##################################################### Day 0 ######################################################################

FS9.D0.g <- subset_samples(FS9.genus, Day == 'D0')
sample_sums(FS9.D0.g)
colnames(otu_table(FS9.D0.g)) #check on all the sample names
FS9.D0.g <- prune_taxa(taxa_sums(FS9.D0.g) > 1, FS9.D0.g)
#if taxa_sums is >1, then it will print that out in FS9.D0.g object and not include anything with <1.
rowSums(FS9.D0.g@otu_table)

#Look at what Set is
unique(sample_data(FS9.D0.g)$Set)

FS9.D0.g.De <- phyloseq_to_deseq2(FS9.D0.g, ~ Set)
FS9.D0.g.De <- DESeq(FS9.D0.g.De, test = "Wald", fitType = "parametric")

######### Day 0 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_INFnm")
#8
sum(meta$Set == "D0_NONINFnm")
#9

resultsNames(FS9.D0.g.De)
#[1] [1] "Intercept"           "Set_D0_NONINFnm_vs_D0_INFnm"
sample_data(FS9.D0.g)$Set <- factor(sample_data(FS9.D0.g)$Set,
                                  levels =c('D0_NONINFnm',"D0_INFnm"))
FS9.D0.g.De <- phyloseq_to_deseq2(FS9.D0.g, ~ Set)
FS9.D0.g.De <- DESeq(FS9.D0.g.De, test = "Wald", fitType = "parametric")
FS9.D0.g.De$Set
resultsNames(FS9.D0.g.De)
#[1] "Intercept"             "Set_D0_INFnm_vs_D0_NONINFnm"
res.D0.g = lfcShrink(FS9.D0.g.De, coef = "Set_D0_INFnm_vs_D0_NONINFnm", type = 'apeglm')
sigtab.D0.g = res.D0.g[which(res.D0.g$padj < .05), ]
sigtab.D0.g = cbind(as(sigtab.D0.g, "data.frame"), as(tax_table(FS9.D0.g)[rownames(sigtab.D0.g), ], "matrix"))
sigtab.D0.g$newp <- format(round(sigtab.D0.g$padj, digits = 3), scientific = TRUE)
sigtab.D0.g$Treatment <- ifelse(sigtab.D0.g$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D0.g)

deseq.D0.g <- 
  ggplot(sigtab.D0.g, aes(x=reorder(rownames(sigtab.D0.g), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.g), y=1, label = paste(Order,Genus, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Order Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day 0')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.D0.g

#Add OTU and comparisons columns
sigtab.D0.g
sigtab.D0.g$OTU <- rownames(sigtab.D0.g)
sigtab.D0.g
sigtab.D0.g$comp <- 'D0_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab.g <- rbind(final.sigtab.g, sigtab.D0.g)

######################################################### Day 4 #########################################################

unique(sample_data(FS9.genus)$Set)

FS9.D4.g <- subset_samples(FS9.genus, Day == 'D4')
sample_sums(FS9.D4.g)
colnames(otu_table(FS9.D4.g)) #check on all the sample names
FS9.D4.g <- prune_taxa(taxa_sums(FS9.D4.g) > 1, FS9.D4.g)
#if taxa_sums is >1, then it will print that out in FS9.D4.g object and not include anything with <1.
rowSums(FS9.D4.g@otu_table)

#Look at what Set is
unique(sample_data(FS9.D4.g)$Set)

FS9.D4.g.De <- phyloseq_to_deseq2(FS9.D4.g, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.D4.g.De <- DESeq(FS9.D4.g.De, test = "Wald", fitType = "parametric")

######### Day 4 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_INFnm")
#6
sum(meta$Set == "D4_NONINFnm")
#5

resultsNames(FS9.D4.g.De)
#[1] "Intercept"             "Set_D4_NONINFnm_vs_D4_INFnm"

#re-level your factor and re-run DESeq2
sample_data(FS9.D4.g)$Set <- factor(sample_data(FS9.D4.g)$Set,
                                  levels =c('D4_NONINFnm',"D4_INFnm"))
FS9.D4.g.De <- phyloseq_to_deseq2(FS9.D4.g, ~ Set)
FS9.D4.g.De <- DESeq(FS9.D4.g.De, test = "Wald", fitType = "parametric")
FS9.D4.g.De$Set
resultsNames(FS9.D4.g.De)
#[1] "Intercept"              "Set_D4_INFnm_vs_D4_NONINFnm"
res.D4.g = lfcShrink(FS9.D4.g.De, coef = "Set_D4_INFnm_vs_D4_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.D4.g = res.D4.g[which(res.D4.g$padj < .05), ]
sigtab.D4.g = cbind(as(sigtab.D4.g, "data.frame"), as(tax_table(FS9.D4.g)[rownames(sigtab.D4.g), ], "matrix"))
sigtab.D4.g$newp <- format(round(sigtab.D4.g$padj, digits = 3), scientific = TRUE)
sigtab.D4.g$Treatment <- ifelse(sigtab.D4.g$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D4.g)
#DataFrame with 0 rows and 7 columns. Will continue to next data set


######################################################### Day 7 #########################################################

unique(sample_data(FS9.genus)$Set)

FS9.D7.g <- subset_samples(FS9.genus, Day == 'D7')
sample_sums(FS9.D7.g)
colnames(otu_table(FS9.D7.g)) #check on all the sample names
FS9.D7.g <- prune_taxa(taxa_sums(FS9.D7.g) > 1, FS9.D7.g)
#if taxa_sums is >1, then it will print that out in FS9.D7.g object and not include anything with <1.
rowSums(FS9.D7.g@otu_table)

#Look at what Set is
unique(sample_data(FS9.D7.g)$Set)

FS9.D7.g.De <- phyloseq_to_deseq2(FS9.D7.g, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.D7.g.De <- DESeq(FS9.D7.g.De, test = "Wald", fitType = "parametric")

######### Day 7 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFnm")
#14
sum(meta$Set == "D7_NONINFnm")
#12

sample_data(FS9.D7.g)$Set <- factor(sample_data(FS9.D7.g)$Set,
                                  levels =c('D7_NONINFnm',"D7_INFnm"))
FS9.D7.g.De <- phyloseq_to_deseq2(FS9.D7.g, ~ Set)
FS9.D7.g.De <- DESeq(FS9.D7.g.De, test = "Wald", fitType = "parametric")
FS9.D7.g.De$Set
resultsNames(FS9.D7.g.De)
#[1] "Intercept"              "Set_D7_INFnm_vs_D7_NONINFnm"
res.D7.g = lfcShrink(FS9.D7.g.De, coef = "Set_D7_INFnm_vs_D7_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.D7.g = res.D7.g[which(res.D7.g$padj < .05), ]
sigtab.D7.g = cbind(as(sigtab.D7.g, "data.frame"), as(tax_table(FS9.D7.g)[rownames(sigtab.D7.g), ], "matrix"))
sigtab.D7.g$newp <- format(round(sigtab.D7.g$padj, digits = 3), scientific = TRUE)
sigtab.D7.g$Treatment <- ifelse(sigtab.D7.g$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D7.g)
#DataFrame with 0 rows and 7 columns

#######################################################################################################

#write csv that includes all differentially abundant OTUs regardless of log2-fold values
write.csv(final.sigtab.g, file= "FS9_FinalDiffAbund_Genus_OutDoubletons_Q1_AllDays.csv")

#######################################################################################################

##################################### Q2: SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (GENUS LEVEL) BETWEEN TREATMENTS ###################################


#Import files
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
#otu_table()   OTU Table:         [ 2039 taxa and 47 samples ]
#sample_data() Sample Data:       [ 47 samples by 5 sample variables ]
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
FS9.genus <- tax_glom(FS9, taxrank = "Genus")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 

##################################### Day 4, 7 ##########################################################################

######################################################### Day 4 #########################################################

sample_data(FS9.genus)

FS9.D4 <- subset_samples(FS9.genus, Day == 'D4')
sample_sums(FS9.D4)
colnames(otu_table(FS9.D4)) #check on all the sample names
FS9.D4 <- prune_taxa(taxa_sums(FS9.D4) > 1, FS9.D4)
#if taxa_sums is >1, then it will print that out in FS9.D4 object and not include anything with <1.
rowSums(FS9.D4@otu_table)

######### 1. Day 4 INFinject vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_INFnm")
#6
sum(meta$Set == "D4_INFinject")
#7

sample_data(FS9.D4)$Set <- factor(sample_data(FS9.D4)$Set,
                                  levels =c('D4_INFnm', "D4_INFinject", "D4_INFfeed"))
FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")
FS9.D4.De$Set
resultsNames(FS9.D4.De)
#[1] "Intercept"        "Set_D4_INFinject_vs_D4_INFnm" "Set_D4_INFfeed_vs_D4_INFnm" 
res.D4.ji = lfcShrink(FS9.D4.De, coef = "Set_D4_INFinject_vs_D4_INFnm", type = 'apeglm')
sigtab.D4.ji = res.D4.ji[which(res.D4.ji$padj < .05), ]
sigtab.D4.ji = cbind(as(sigtab.D4.ji, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.ji), ], "matrix"))
sigtab.D4.ji$newp <- format(round(sigtab.D4.ji$padj, digits = 3), scientific = TRUE)
sigtab.D4.ji$Treatment <- ifelse(sigtab.D4.ji$log2FoldChange >=0, "INFinject", "INFnm")
head(sigtab.D4.ji) #1 OTU

deseq.D4.ji <- 
  ggplot(sigtab.D4.ji, aes(x=reorder(rownames(sigtab.D4.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.ji), y=3, label = paste(Order,Genus, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Order Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Genus\n in INFinject Group Relative to INFnm\n in Fecal Microbiota on Day 4')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066'))
deseq.D4.ji

#Add OTU and comparisons columns
sigtab.D4.ji
sigtab.D4.ji$OTU <- rownames(sigtab.D4.ji)
sigtab.D4.ji
sigtab.D4.ji$comp <- 'D4_INFinject_vs_INFnm'

#Create final significant comparisons table
final.sigtab <- sigtab.D4.ji


######### 2. Day 4 INFfeed vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_INFnm")
#6
sum(meta$Set == "D4_INFfeed")
#4

#Extract results from a DESeq analysis, organize table
FS9.D4.De$Set
resultsNames(FS9.D4.De)
#[1] "Intercept"         "Set_D4_INFinject_vs_D4_INFnm" "Set_D4_INFfeed_vs_D4_INFnm" 
res.D4.oi = lfcShrink(FS9.D4.De, coef = "Set_D4_INFfeed_vs_D4_INFnm", type = 'apeglm')
sigtab.D4.oi = res.D4.oi[which(res.D4.oi$padj < .05), ]
sigtab.D4.oi = cbind(as(sigtab.D4.oi, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.oi), ], "matrix"))
format(sigtab.D4.oi$padj, scientific = TRUE)
sigtab.D4.oi$newp <- format(round(sigtab.D4.oi$padj, digits = 3), scientific = TRUE)
sigtab.D4.oi$Treatment <- ifelse(sigtab.D4.oi$log2FoldChange >=0, "INFfeed", "INFnm")
head(sigtab.D4.oi)  #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


######### 3. Day 4 INFfeed vs INFinject ###################

#INFnm = I
#INFinject= J
#INFfeed = O

#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_INFinject")
#7
sum(meta$Set == "D4_INFfeed")
#4

#Extract results from a DESeq analysis, organize table
resultsNames(FS9.D4.De)
#[1] "Intercept"           "Set_D4_INFinject_vs_D4_INFnm" "Set_D4_INFfeed_vs_D4_INFnm"  
sample_data(FS9.D4)$Set <- factor(sample_data(FS9.D4)$Set,
                                  levels =c("D4_INFinject", "D4_INFfeed",
                                            'D4_INFnm'))
FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D4.De)
#[1] "Intercept"            "Set_D4_INFfeed_vs_D4_INFinject" "Set_D4_INFnm_vs_D4_INFinject" 
res.D4.oj = lfcShrink(FS9.D4.De, coef = "Set_D4_INFfeed_vs_D4_INFinject", type = 'apeglm')
sigtab.D4.oj = res.D4.oi[which(res.D4.oi$padj < .05), ]
sigtab.D4.oj = cbind(as(sigtab.D4.oj, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.oj), ], "matrix"))
format(sigtab.D4.oj$padj, scientific = TRUE)
sigtab.D4.oj$newp <- format(round(sigtab.D4.oj$padj, digits = 3), scientific = TRUE)
sigtab.D4.oj$Treatment <- ifelse(sigtab.D4.oj$log2FoldChange >=0, "INFfeed", "INFinject")
head(sigtab.D4.oj) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison



##################################################### Day 7 ######################################################################

sample_data(FS9.genus)
FS9.D7 <- subset_samples(FS9.genus, Day == 'D7')
sample_sums(FS9.D7)
colnames(otu_table(FS9.D7)) #check on all the sample names
FS9.D7 <- prune_taxa(taxa_sums(FS9.D7) > 1, FS9.D7)
#if taxa_sums is >1, then it will print that out in FS9.D7 object and not include anything with <1.
rowSums(FS9.D7@otu_table)

#Look at what Set is
sample_data(FS9.D7)

######### 1. Day 7 INFinject vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFnm")
#14
sum(meta$Set == "D7_INFinject")
#19

sample_data(FS9.D7)$Set <- factor(sample_data(FS9.D7)$Set,
                                  levels =c("D7_INFinject", 'D7_INFnm', "D7_INFfeed"))
FS9.D7.De <- phyloseq_to_deseq2(FS9.D7, ~ Set)
FS9.D7.De <- DESeq(FS9.D7.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.De)
#[1]  "Intercept"        "Set_D7_INFnm_vs_D7_INFinject"   "Set_D7_INFfeed_vs_D7_INFinject" 
res.D7.ji = lfcShrink(FS9.D7.De, coef = "Set_D7_INFnm_vs_D7_INFinject", type = 'apeglm')
sigtab.D7.ji = res.D7.ji[which(res.D7.ji$padj < .05), ]
sigtab.D7.ji = cbind(as(sigtab.D7.ji, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7.ji), ], "matrix"))
sigtab.D7.ji$newp <- format(round(sigtab.D7.ji$padj, digits = 3), scientific = TRUE)
sigtab.D7.ji$Treatment <- ifelse(sigtab.D7.ji$log2FoldChange >=0, "Depleted in INFinject", "Enriched in INFinject")
head(sigtab.D7.ji)

deseq.D7.ji <- 
  ggplot(sigtab.D7.ji, aes(x=reorder(rownames(sigtab.D7.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.ji), y=3, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFinject Group Relative to INFnm\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066'))
deseq.D7.ji

#Add OTU and comparisons columns
sigtab.D7.ji
sigtab.D7.ji$OTU <- rownames(sigtab.D7.ji)
sigtab.D7.ji
sigtab.D7.ji$comp <- 'D7_INFnm_vs_INFinject'

#Create final significant comparisons table
final.sigtab <- sigtab.D7.ji



######### 2. Day 7 INFfeed vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFnm")
#14
sum(meta$Set == "D7_INFfeed")
#13

#Extract results from a DESeq analysis, organize table
sample_data(FS9.D7)$Set <- factor(sample_data(FS9.D7)$Set,
                                  levels =c("D7_INFfeed", 'D7_INFnm', 'D7_INFinject'))
FS9.D7.De <- phyloseq_to_deseq2(FS9.D7, ~ Set)
FS9.D7.De <- DESeq(FS9.D7.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.De)
#[1] "Intercept"               "Set_D7_INFnm_vs_D7_INFfeed"     "Set_D7_INFinject_vs_D7_INFfeed"   
res.D7.oi = lfcShrink(FS9.D7.De, coef = "Set_D7_INFnm_vs_D7_INFfeed", type = 'apeglm')
sigtab.D7.oi = res.D7.oi[which(res.D7.oi$padj < .05), ]
sigtab.D7.oi = cbind(as(sigtab.D7.oi, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7.oi), ], "matrix"))
format(sigtab.D7.oi$padj, scientific = TRUE)
sigtab.D7.oi$newp <- format(round(sigtab.D7.oi$padj, digits = 3), scientific = TRUE)
sigtab.D7.oi$Treatment <- ifelse(sigtab.D7.oi$log2FoldChange >=0, "Depleted in INFfeed", "Enriched in INFfeed")

#Summarize sigtab.D7.oi
sum.sigtab.D7.oi <- summary(sigtab.D7.oi)
sum.sigtab.D7.oi

#ggplot
deseq.D7.oi <- ggplot(sigtab.D7.oi, aes(x=reorder(rownames(sigtab.D7.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.oi), y=3, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFfeed Group Relative to INFnm\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', INFnm='#CC0066'))
deseq.D7.oi

#Add OTU and comparisons columns
sigtab.D7.oi
sigtab.D7.oi$OTU <- rownames(sigtab.D7.oi)
sigtab.D7.oi
sigtab.D7.oi$comp <- 'D7_INFfeed_vs_INFnm'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.D7.oi, final.sigtab)