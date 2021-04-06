#####################################################################################################
#FS9 DESeq2 - Infected: nm vs feed vs inject, days 7, 11, 14 - order level
#By Mou, KT

#Purpose: This code uses DESeq2 package to identify fecal microbial genera that were differentially 
#abundant between the three groups on days 7, 11, 14

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


##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (ORDER LEVEL) BETWEEN TREATMENTS #############

##################################################### Day 7 ######################################################################

sample_data(FS9.order)
FS9.D7 <- subset_samples(FS9.order, Day == 'D7')
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
#8
sum(meta$Set == "D7_INFinject")
#13

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
#DataFrame with 0 rows and 7 columns

######### 2. Day 7 INFfeed vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFnm")
#8
sum(meta$Set == "D7_INFfeed")
#5

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
sigtab.D7.oi
#DataFrame with 0 rows and 7 columns

##################################################### Day 11 ######################################################################

sample_data(FS9.order)

FS9.D11 <- subset_samples(FS9.order, Day == 'D11')
sample_sums(FS9.D11)
colnames(otu_table(FS9.D11)) #check on all the sample names
FS9.D11 <- prune_taxa(taxa_sums(FS9.D11) > 1, FS9.D11)
#if taxa_sums is >1, then it will print that out in FS9.D11 object and not include anything with <1.
rowSums(FS9.D11@otu_table)

######### 1. Day 11 INFinject vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D11_INFnm")
#6
sum(meta$Set == "D11_INFinject")
#7

sample_data(FS9.D11)$Set <- factor(sample_data(FS9.D11)$Set,
                                  levels =c("D11_INFinject", 'D11_INFnm', "D11_INFfeed"))
FS9.D11.De <- phyloseq_to_deseq2(FS9.D11, ~ Set)
FS9.D11.De <- DESeq(FS9.D11.De, test = "Wald", fitType = "parametric")
FS9.D11.De$Set
resultsNames(FS9.D11.De)
#[1] "Intercept"        "Set_D11_INFinject_vs_D11_INFnm" "Set_D11_INFfeed_vs_D11_INFnm" 
res.D11.ji = lfcShrink(FS9.D11.De, coef = "Set_D11_INFnm_vs_D11_INFinject", type = 'apeglm')
sigtab.D11.ji = res.D11.ji[which(res.D11.ji$padj < .05), ]
sigtab.D11.ji = cbind(as(sigtab.D11.ji, "data.frame"), as(tax_table(FS9.D11)[rownames(sigtab.D11.ji), ], "matrix"))
sigtab.D11.ji$newp <- format(round(sigtab.D11.ji$padj, digits = 3), scientific = TRUE)
sigtab.D11.ji$Treatment <- ifelse(sigtab.D11.ji$log2FoldChange >=0, "Depleted in INFinject", "Enriched in INFinject")
head(sigtab.D11.ji) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


######### 2. Day 11 INFfeed vs INFnm ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D11_INFnm")
#6
sum(meta$Set == "D11_INFfeed")
#4

#Extract results from a DESeq analysis, organize table
FS9.D11.De$Set
sample_data(FS9.D11)$Set <- factor(sample_data(FS9.D11)$Set,
                                   levels =c("D11_INFfeed", 'D11_INFnm', "D11_INFinject"))
FS9.D11.De <- phyloseq_to_deseq2(FS9.D11, ~ Set)
FS9.D11.De <- DESeq(FS9.D11.De, test = "Wald", fitType = "parametric")
FS9.D11.De$Set
resultsNames(FS9.D11.De)
#[1] "Intercept"         "Set_D11_INFnm_vs_D11_INFfeed"     "Set_D11_INFinject_vs_D11_INFfeed"
res.D11.oi = lfcShrink(FS9.D11.De, coef = "Set_D11_INFnm_vs_D11_INFfeed", type = 'apeglm')
sigtab.D11.oi = res.D11.oi[which(res.D11.oi$padj < .05), ]
sigtab.D11.oi = cbind(as(sigtab.D11.oi, "data.frame"), as(tax_table(FS9.D11)[rownames(sigtab.D11.oi), ], "matrix"))
format(sigtab.D11.oi$padj, scientific = TRUE)
sigtab.D11.oi$newp <- format(round(sigtab.D11.oi$padj, digits = 3), scientific = TRUE)
sigtab.D11.oi$Treatment <- ifelse(sigtab.D11.oi$log2FoldChange >=0, "Depleted in INFfeed", "Enriched in INFfeed")
head(sigtab.D11.oi)  #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison

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



################### Create plot containing all significant orders with log2-fold changes > 0.25 for DeSeq2 results

#ggplot
deseqfinalplot <- ggplot(final.sigtab, aes(x=reorder(rownames(final.sigtab), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(final.sigtab), y=3, label = paste(Order)), size=5, fontface = 'italic')+ labs(x="Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12)) +
  coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c('Depleted in INFfeed'='dodgerblue1', 'Depleted in INFinject'='green3')) +
  theme_bw()
#light green for depleted in INFinject
#light blue for INFfeed
#color name options for ggplot: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
deseqfinalplot

ggsave("Q2_INFnm_INFinject_INFfeed_Order_D14_DeSeq2_DepletedPlots.tiff", plot=deseqfinalplot, width = 8, height = 4, dpi = 500, units =c("in"))

#write csv
write.csv(final.sigtab, file= "FS9_FinalDiffAbund_Order_OutDoubletons_Q2_D47D11D14.csv")
