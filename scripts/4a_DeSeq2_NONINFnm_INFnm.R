#####################################################################################################
#FS9 DESeq2 - Noninfected (NONINFnm) vs Infected (INFnm) only on all days, or days 7, 11, and 14
#By Mou, KT

#Purpose: This code uses DESeq2 package to identify fecal microbial taxa that were differentially 
#abundant between the two groups NONINFnm and INFnm on all days and only days 7, 11, and 14 at the genus and order level

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

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (ORDER LEVEL) BETWEEN TREATMENTS #############

##################################### Days 7, 11, and 14 ##########################################################################

otu <- import_mothur(mothur_shared_file = './data/stability.outdoubletons.abund.opti_mcc.shared') #use unrarified data
taxo <- import_mothur(mothur_constaxonomy_file = './data/stability.outdoubletons.abund.opti_mcc.0.03.cons.taxonomy')
meta <- read.table(file = './data/FS9_metadata_NONINFnm_vs_INFnm_alldays.csv', sep = ',', header = TRUE)

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
#otu_table()   OTU Table:         [ 2039 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 5 sample variables ]
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

######### Day 7 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFnm")
#8
sum(meta$Set == "D7_NONINFnm")
#9

sample_data(FS9.D7)$Set <- factor(sample_data(FS9.D7)$Set,
                                     levels =c('D7_NONINFnm',"D7_INFnm"))
FS9.D7.De <- phyloseq_to_deseq2(FS9.D7, ~ Set)
FS9.D7.De <- DESeq(FS9.D7.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.De)
#[1] "Intercept"                 "Set_D7_INFnm_vs_D7_NONINFnm"
res.D7 = lfcShrink(FS9.D7.De, coef = "Set_D7_INFnm_vs_D7_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.D7 = res.D7[which(res.D7$padj < .05), ]
sigtab.D7 = cbind(as(sigtab.D7, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7), ], "matrix"))
sigtab.D7$newp <- format(round(sigtab.D7$padj, digits = 3), scientific = TRUE)
sigtab.D7$Treatment <- ifelse(sigtab.D7$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D7)
#DataFrame with 0 rows and 7 columns. Skip to next data set


##################################################### Day 11 ######################################################################

FS9.D11 <- subset_samples(FS9.order, Day == 'D11')
sample_sums(FS9.D11)
colnames(otu_table(FS9.D11)) #check on all the sample names
FS9.D11 <- prune_taxa(taxa_sums(FS9.D11) > 1, FS9.D11)
#if taxa_sums is >1, then it will print that out in FS9.D11 object and not include anything with <1.
rowSums(FS9.D11@otu_table)

#Look at what Set is
unique(sample_data(FS9.D11)$Set)

######### Day 11 INFnm vs NONINFnm ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D11_INFnm")
#6
sum(meta$Set == "D11_NONINFnm")
#5

sample_data(FS9.D11)$Set <- factor(sample_data(FS9.D11)$Set,
                                  levels =c('D11_NONINFnm',"D11_INFnm"))
FS9.D11.De <- phyloseq_to_deseq2(FS9.D11, ~ Set)
FS9.D11.De <- DESeq(FS9.D11.De, test = "Wald", fitType = "parametric")
FS9.D11.De$Set
resultsNames(FS9.D11.De)
#[1] "Intercept"           "Set_D11_INFnm_vs_D11_NONINFnm"
res.D11 = lfcShrink(FS9.D11.De, coef = "Set_D11_INFnm_vs_D11_NONINFnm", type = 'apeglm')
sigtab.D11 = res.D11[which(res.D11$padj < .05), ]
sigtab.D11 = cbind(as(sigtab.D11, "data.frame"), as(tax_table(FS9.D11)[rownames(sigtab.D11), ], "matrix"))
sigtab.D11$newp <- format(round(sigtab.D11$padj, digits = 3), scientific = TRUE)
sigtab.D11$Treatment <- ifelse(sigtab.D11$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D11)
#DataFrame with 0 rows and 7 columns. Skip to next data set

##################################################### Day 14 ######################################################################

FS9.D14 <- subset_samples(FS9.order, Day == 'D14')
sample_sums(FS9.D14)
colnames(otu_table(FS9.D14)) #check on all the sample names
FS9.D14 <- prune_taxa(taxa_sums(FS9.D14) > 1, FS9.D14)
#if taxa_sums is >1, then it will print that out in FS9.D14 object and not include anything with <1.
rowSums(FS9.D14@otu_table)

#Look at what Set is
unique(sample_data(FS9.D14)$Set)

######### Day 14 INFnm vs NONINFnm ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D14_INFnm")
#14
sum(meta$Set == "D14_NONINFnm")
#12

sample_data(FS9.D14)$Set <- factor(sample_data(FS9.D14)$Set,
                                   levels =c('D14_NONINFnm',"D14_INFnm"))
FS9.D14.De <- phyloseq_to_deseq2(FS9.D14, ~ Set)
FS9.D14.De <- DESeq(FS9.D14.De, test = "Wald", fitType = "parametric")
FS9.D14.De$Set
resultsNames(FS9.D14.De)
#[1] "Intercept"           "Set_D14_INFnm_vs_D14_NONINFnm"
res.D14 = lfcShrink(FS9.D14.De, coef = "Set_D14_INFnm_vs_D14_NONINFnm", type = 'apeglm')
sigtab.D14 = res.D14[which(res.D14$padj < .05), ]
sigtab.D14 = cbind(as(sigtab.D14, "data.frame"), as(tax_table(FS9.D14)[rownames(sigtab.D14), ], "matrix"))
sigtab.D14$newp <- format(round(sigtab.D14$padj, digits = 3), scientific = TRUE)
sigtab.D14$Treatment <- ifelse(sigtab.D14$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D14)
#baseMean log2FoldChange       lfcSE       pvalue       padj  Kingdom          Phylum            Class              Order Family
#Otu0112 43.29473   8.157507e-07 0.001442695 2.880172e-05 0.00100806 Bacteria Verrucomicrobia Verrucomicrobiae Verrucomicrobiales   <NA>
#  Genus  newp Treatment
#Otu0112  <NA> 1e-03     INFnm

sigtab.D14
sigtab.D14$OTU <- rownames(sigtab.D14)
sigtab.D14
sigtab.D14$comp <- 'D14_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- sigtab.D14

#Log2-fold change is too small. Not worth mentioning.

#######################################################################################################

#write csv that includes all differentially abundant OTUs regardless of log2-fold values
write.csv(final.sigtab, file= "FS9_FinalDiffAbund_Order_OutDoubletons_Q1_D7D11D14.csv")

#######################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (GENUS LEVEL) BETWEEN TREATMENTS #########

##################################### Day 7, 11, 14 ##########################################################################

##################################################### Day 7 ######################################################################

FS9.D7.g <- subset_samples(FS9.genus, Day == 'D7')
sample_sums(FS9.D7.g)
colnames(otu_table(FS9.D7.g)) #check on all the sample names
FS9.D7.g <- prune_taxa(taxa_sums(FS9.D7.g) > 1, FS9.D7.g)
#if taxa_sums is >1, then it will print that out in FS9.D7.g object and not include anything with <1.
rowSums(FS9.D7.g@otu_table)

#Look at what Set is
unique(sample_data(FS9.D7.g)$Set)

FS9.D7.g.De <- phyloseq_to_deseq2(FS9.D7.g, ~ Set)
FS9.D7.g.De <- DESeq(FS9.D7.g.De, test = "Wald", fitType = "parametric")

######### Day 7 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFnm")
#8
sum(meta$Set == "D7_NONINFnm")
#9

resultsNames(FS9.D7.g.De)
#[1] [1] "Intercept"           "Set_D7_NONINFnm_vs_D7_INFnm"
sample_data(FS9.D7.g)$Set <- factor(sample_data(FS9.D7.g)$Set,
                                    levels =c('D7_NONINFnm',"D7_INFnm"))
FS9.D7.g.De <- phyloseq_to_deseq2(FS9.D7.g, ~ Set)
FS9.D7.g.De <- DESeq(FS9.D7.g.De, test = "Wald", fitType = "parametric")
FS9.D7.g.De$Set
resultsNames(FS9.D7.g.De)
#[1] "Intercept"             "Set_D7_INFnm_vs_D7_NONINFnm"
res.D7.g = lfcShrink(FS9.D7.g.De, coef = "Set_D7_INFnm_vs_D7_NONINFnm", type = 'apeglm')
sigtab.D7.g = res.D7.g[which(res.D7.g$padj < .05), ]
sigtab.D7.g = cbind(as(sigtab.D7.g, "data.frame"), as(tax_table(FS9.D7.g)[rownames(sigtab.D7.g), ], "matrix"))
sigtab.D7.g$newp <- format(round(sigtab.D7.g$padj, digits = 3), scientific = TRUE)
sigtab.D7.g$Treatment <- ifelse(sigtab.D7.g$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D7.g)

deseq.D7.g <- 
  ggplot(sigtab.D7.g, aes(x=reorder(rownames(sigtab.D7.g), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.g), y=0, label = paste(Genus, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Order Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Genus\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#619CFF', NONINFnm="#C77CFF"))
deseq.D7.g

#Add OTU and comparisons columns
sigtab.D7.g
sigtab.D7.g$OTU <- rownames(sigtab.D7.g)
sigtab.D7.g
sigtab.D7.g$comp <- 'D7_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab.g <- sigtab.D7.g

######################################################### Day 11 #########################################################

unique(sample_data(FS9.genus)$Set)

FS9.D11.g <- subset_samples(FS9.genus, Day == 'D11')
sample_sums(FS9.D11.g)
colnames(otu_table(FS9.D11.g)) #check on all the sample names
FS9.D11.g <- prune_taxa(taxa_sums(FS9.D11.g) > 1, FS9.D11.g)
#if taxa_sums is >1, then it will print that out in FS9.D11.g object and not include anything with <1.
rowSums(FS9.D11.g@otu_table)

#Look at what Set is
unique(sample_data(FS9.D11.g)$Set)

FS9.D11.g.De <- phyloseq_to_deseq2(FS9.D11.g, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.D11.g.De <- DESeq(FS9.D11.g.De, test = "Wald", fitType = "parametric")

######### Day 11 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D11_INFnm")
#6
sum(meta$Set == "D11_NONINFnm")
#5

resultsNames(FS9.D11.g.De)
#[1] "Intercept"             "Set_D11_NONINFnm_vs_D11_INFnm"

#re-level your factor and re-run DESeq2
sample_data(FS9.D11.g)$Set <- factor(sample_data(FS9.D11.g)$Set,
                                    levels =c('D11_NONINFnm',"D11_INFnm"))
FS9.D11.g.De <- phyloseq_to_deseq2(FS9.D11.g, ~ Set)
FS9.D11.g.De <- DESeq(FS9.D11.g.De, test = "Wald", fitType = "parametric")
FS9.D11.g.De$Set
resultsNames(FS9.D11.g.De)
#[1] "Intercept"              "Set_D11_INFnm_vs_D11_NONINFnm"
res.D11.g = lfcShrink(FS9.D11.g.De, coef = "Set_D11_INFnm_vs_D11_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.D11.g = res.D11.g[which(res.D11.g$padj < .05), ]
sigtab.D11.g = cbind(as(sigtab.D11.g, "data.frame"), as(tax_table(FS9.D11.g)[rownames(sigtab.D11.g), ], "matrix"))
sigtab.D11.g$newp <- format(round(sigtab.D11.g$padj, digits = 3), scientific = TRUE)
sigtab.D11.g$Treatment <- ifelse(sigtab.D11.g$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D11.g)
#DataFrame with 0 rows and 7 columns. Will continue to next data set


######################################################### Day 14 #########################################################

unique(sample_data(FS9.genus)$Set)

FS9.D14.g <- subset_samples(FS9.genus, Day == 'D14')
sample_sums(FS9.D14.g)
colnames(otu_table(FS9.D14.g)) #check on all the sample names
FS9.D14.g <- prune_taxa(taxa_sums(FS9.D14.g) > 1, FS9.D14.g)
#if taxa_sums is >1, then it will print that out in FS9.D14.g object and not include anything with <1.
rowSums(FS9.D14.g@otu_table)

#Look at what Set is
unique(sample_data(FS9.D14.g)$Set)

FS9.D14.g.De <- phyloseq_to_deseq2(FS9.D14.g, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with
# Set combined day and treatment

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
FS9.D14.g.De <- DESeq(FS9.D14.g.De, test = "Wald", fitType = "parametric")

######### Day 14 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D14_INFnm")
#14
sum(meta$Set == "D14_NONINFnm")
#12

sample_data(FS9.D14.g)$Set <- factor(sample_data(FS9.D14.g)$Set,
                                    levels =c('D14_NONINFnm',"D14_INFnm"))
FS9.D14.g.De <- phyloseq_to_deseq2(FS9.D14.g, ~ Set)
FS9.D14.g.De <- DESeq(FS9.D14.g.De, test = "Wald", fitType = "parametric")
FS9.D14.g.De$Set
resultsNames(FS9.D14.g.De)
#[1] "Intercept"              "Set_D14_INFnm_vs_D14_NONINFnm"
res.D14.g = lfcShrink(FS9.D14.g.De, coef = "Set_D14_INFnm_vs_D14_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.D14.g = res.D14.g[which(res.D14.g$padj < .05), ]
sigtab.D14.g = cbind(as(sigtab.D14.g, "data.frame"), as(tax_table(FS9.D14.g)[rownames(sigtab.D14.g), ], "matrix"))
sigtab.D14.g$newp <- format(round(sigtab.D14.g$padj, digits = 3), scientific = TRUE)
sigtab.D14.g$Treatment <- ifelse(sigtab.D14.g$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D14.g)
#DataFrame with 0 rows and 7 columns

#######################################################################################################

#write csv that includes all differentially abundant OTUs regardless of log2-fold values
write.csv(final.sigtab.g, file= "FS9_FinalDiffAbund_Genus_OutDoubletons_Q1_D7D11D14.csv")

#######################################################################################################