#####################################################################################################
#FS9 DESeq2 - Infected: nm vs feed vs inject, days 4 and 7
#Kathy Mou

#Purpose: This code uses DESeq2 package to identify fecal microbial genera that were differentially 
#abundant between the three groups on days 4 and 7

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
FS9.order <- tax_glom(FS9, taxrank = "Order")
FS9.phylum <- tax_glom(FS9, taxrank = "Phylum")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 


##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (ORDER LEVEL) BETWEEN TREATMENTS ###################################

##################################################### Day 4 ######################################################################


sample_data(FS9.order)

FS9.D4 <- subset_samples(FS9.order, Day == 'D4')
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
head(sigtab.D4.ji) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


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
sigtab.D7.ji$Treatment <- ifelse(sigtab.D7.ji$log2FoldChange >=0, "INFnm", "INFinject")
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
sigtab.D7.oi$Treatment <- ifelse(sigtab.D7.oi$log2FoldChange >=0, "INFnm", "INFfeed")

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


######### 3. Day 7 INFfeed vs INFinject ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFinject")
#19
sum(meta$Set == "D7_INFfeed")
#13

#Extract results from a DESeq analysis, organize table
sample_data(FS9.D7)$Set <- factor(sample_data(FS9.D7)$Set,
                                     levels =c("D7_INFfeed", "D7_INFinject", 
                                               'D7_INFnm'))
FS9.D7.De <- phyloseq_to_deseq2(FS9.D7, ~ Set)
FS9.D7.De <- DESeq(FS9.D7.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.De)
#[1] "Intercept"             "Set_D7_INFinject_vs_D7_INFfeed" "Set_D7_INFnm_vs_D7_INFfeed" 
res.D7.oj = lfcShrink(FS9.D7.De, coef = "Set_D7_INFinject_vs_D7_INFfeed", type = 'apeglm')
sigtab.D7.oj = res.D7.oi[which(res.D7.oi$padj < .05), ]
sigtab.D7.oj = cbind(as(sigtab.D7.oj, "data.frame"), as(tax_table(FS9.D7)[rownames(sigtab.D7.oj), ], "matrix"))
format(sigtab.D7.oj$padj, scientific = TRUE)
sigtab.D7.oj$newp <- format(round(sigtab.D7.oj$padj, digits = 3), scientific = TRUE)
sigtab.D7.oj$Treatment <- ifelse(sigtab.D7.oj$log2FoldChange >=0, "INFinject", "INFfeed")

#Summarize sigtab.D7.oj
sum.sigtab.D7.oj <- summary(sigtab.D7.oj)
sum.sigtab.D7.oj

#ggplot
deseq.D7.oj <- ggplot(sigtab.D7.oj, aes(x=reorder(rownames(sigtab.D7.oj), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.oj), y=3, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Order\n in INFfeed Group Relative to INFinject\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', INFinject='#E69F00'))
deseq.D7.oj

#Add OTU and comparisons columns
sigtab.D7.oj
sigtab.D7.oj$OTU <- rownames(sigtab.D7.oj)
sigtab.D7.oj
sigtab.D7.oj$comp <- 'D7_INFfeed_vs_INFinject'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.D7.oj, final.sigtab)






#######################################################################################################

#write csv
write.csv(final.sigtab, file= "FS9_FinalDiffAbund_Order_OutDoubletons_Q2_D4D7.csv")

#######################################################################################################






##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (PHYLUM LEVEL) BETWEEN TREATMENTS ###################################

##################################################### Day 4 ######################################################################

FS9.D4.p <- subset_samples(FS9.phylum, Day == 'D4')
sample_sums(FS9.D4.p)
colnames(otu_table(FS9.D4.p)) #check on all the sample names
FS9.D4.p <- prune_taxa(taxa_sums(FS9.D4.p) > 1, FS9.D4.p)
#if taxa_sums is >1, then it will print that out in FS9.D4.p object and not include anything with <1.
rowSums(FS9.D4.p@otu_table)

######### 1. Day 4 INFinject vs INFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_INFnm")
#6
sum(meta$Set == "D4_INFinject")
#7

sample_data(FS9.D4.p)$Set <- factor(sample_data(FS9.D4.p)$Set,
                                  levels =c('D4_INFnm', "D4_INFinject", "D4_INFfeed"))
FS9.D4.p.De <- phyloseq_to_deseq2(FS9.D4.p, ~ Set)
FS9.D4.p.De <- DESeq(FS9.D4.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D4.p.De)
#[1] "Intercept"          "Set_D4_INFinject_vs_D4_INFnm" "Set_D4_INFfeed_vs_D4_INFnm"  
res.D4.p.ji = lfcShrink(FS9.D4.p.De, coef = "Set_D4_INFinject_vs_D4_INFnm", type = 'apeglm')
sigtab.D4.p.ji = res.D4.p.ji[which(res.D4.p.ji$padj < .05), ]
sigtab.D4.p.ji = cbind(as(sigtab.D4.p.ji, "data.frame"), as(tax_table(FS9.D4.p)[rownames(sigtab.D4.p.ji), ], "matrix"))
sigtab.D4.p.ji$newp <- format(round(sigtab.D4.p.ji$padj, digits = 3), scientific = TRUE)
sigtab.D4.p.ji$Treatment <- ifelse(sigtab.D4.p.ji$log2FoldChange >=0, "INFinject", "INFnm")
head(sigtab.D4.p.ji) #DataFrame with 0 rows and 7 columns, meaning there were no phyla that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


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
FS9.D4.p.De$Set
resultsNames(FS9.D4.p.De)
#[1] "Intercept"           "Set_D4_INFinject_vs_D4_INFnm" "Set_D4_INFfeed_vs_D4_INFnm"  
res.D4.p.oi = lfcShrink(FS9.D4.p.De, coef = "Set_D4_INFfeed_vs_D4_INFnm", type = 'apeglm')
sigtab.D4.p.oi = res.D4.p.oi[which(res.D4.p.oi$padj < .05), ]
sigtab.D4.p.oi = cbind(as(sigtab.D4.p.oi, "data.frame"), as(tax_table(FS9.D4.p)[rownames(sigtab.D4.p.oi), ], "matrix"))
format(sigtab.D4.p.oi$padj, scientific = TRUE)
sigtab.D4.p.oi$newp <- format(round(sigtab.D4.p.oi$padj, digits = 3), scientific = TRUE)
sigtab.D4.p.oi$Treatment <- ifelse(sigtab.D4.p.oi$log2FoldChange >=0, "INFfeed", "INFnm")
head(sigtab.D4.p.oi)  #DataFrame with 0 rows and 7 columns, meaning there were no phyla that were significantly different
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
sample_data(FS9.D4.p)$Set <- factor(sample_data(FS9.D4.p)$Set,
                                  levels =c("D4_INFinject", "D4_INFfeed",
                                            'D4_INFnm'))
FS9.D4.p.De <- phyloseq_to_deseq2(FS9.D4.p, ~ Set)
FS9.D4.p.De <- DESeq(FS9.D4.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D4.p.De)
#[1] "Intercept"                "Set_D4_INFfeed_vs_D4_INFinject" "Set_D4_INFnm_vs_D4_INFinject"
res.D4.p.oj = lfcShrink(FS9.D4.p.De, coef = "Set_D4_INFfeed_vs_D4_INFinject", type = 'apeglm')
sigtab.D4.p.oj = res.D4.p.oi[which(res.D4.p.oi$padj < .05), ]
sigtab.D4.p.oj = cbind(as(sigtab.D4.p.oj, "data.frame"), as(tax_table(FS9.D4.p)[rownames(sigtab.D4.p.oj), ], "matrix"))
format(sigtab.D4.p.oj$padj, scientific = TRUE)
sigtab.D4.p.oj$newp <- format(round(sigtab.D4.p.oj$padj, digits = 3), scientific = TRUE)
sigtab.D4.p.oj$Treatment <- ifelse(sigtab.D4.p.oj$log2FoldChange >=0, "INFfeed", "INFinject")
head(sigtab.D4.p.oj) #DataFrame with 0 rows and 7 columns, meaning there were no phyla that were significantly different
#in abundance between the two groups, so I will skip to the next comparison



##################################################### Day 7 ######################################################################

sample_data(FS9.phylum)
FS9.D7.p <- subset_samples(FS9.phylum, Day == 'D7')
sample_sums(FS9.D7.p)
colnames(otu_table(FS9.D7.p)) #check on all the sample names
FS9.D7.p <- prune_taxa(taxa_sums(FS9.D7.p) > 1, FS9.D7.p)
#if taxa_sums is >1, then it will print that out in FS9.D7.p object and not include anything with <1.
rowSums(FS9.D7.p@otu_table)

#Look at what Set is
sample_data(FS9.D7.p)


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

sample_data(FS9.D7.p)$Set <- factor(sample_data(FS9.D7.p)$Set,
                                  levels =c('D7_INFnm', "D7_INFinject", "D7_INFfeed"))
FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.p.De)
#[1]  "Intercept"             "Set_D7_INFinject_vs_D7_INFnm" "Set_D7_INFfeed_vs_D7_INFnm"    
res.D7.p.ji = lfcShrink(FS9.D7.p.De, coef = "Set_D7_INFinject_vs_D7_INFnm", type = 'apeglm')
sigtab.D7.p.ji = res.D7.p.ji[which(res.D7.p.ji$padj < .05), ]
sigtab.D7.p.ji = cbind(as(sigtab.D7.p.ji, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p.ji), ], "matrix"))
sigtab.D7.p.ji$newp <- format(round(sigtab.D7.p.ji$padj, digits = 3), scientific = TRUE)
sigtab.D7.p.ji$Treatment <- ifelse(sigtab.D7.p.ji$log2FoldChange >=0, "INFinject", "INFnm")
head(sigtab.D7.p.ji)

deseq.D7.p.ji <- 
  ggplot(sigtab.D7.p.ji, aes(x=reorder(rownames(sigtab.D7.p.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.ji), y=-1e-06, label = paste(Phylum, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phylum\n in INFinject Group Relative to INFnm\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066'))
deseq.D7.p.ji

#Add OTU and comparisons columns
sigtab.D7.p.ji
sigtab.D7.p.ji$OTU <- rownames(sigtab.D7.p.ji)
sigtab.D7.p.ji
sigtab.D7.p.ji$comp <- 'D7_INFinject_vs_INFnm'

#Create final significant comparisons table
final.sigtab.phylum <- sigtab.D7.p.ji



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
resultsNames(FS9.D7.p.De)
#[1] "Intercept"      "Set_D7_INFinject_vs_D7_INFnm" "Set_D7_INFfeed_vs_D7_INFnm"  
res.D7.p.oi = lfcShrink(FS9.D7.p.De, coef = "Set_D7_INFfeed_vs_D7_INFnm", type = 'apeglm')
sigtab.D7.p.oi = res.D7.p.oi[which(res.D7.p.oi$padj < .05), ]
sigtab.D7.p.oi = cbind(as(sigtab.D7.p.oi, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p.oi), ], "matrix"))
format(sigtab.D7.p.oi$padj, scientific = TRUE)
sigtab.D7.p.oi$newp <- format(round(sigtab.D7.p.oi$padj, digits = 3), scientific = TRUE)
sigtab.D7.p.oi$Treatment <- ifelse(sigtab.D7.p.oi$log2FoldChange >=0, "INFfeed", "INFnm")
head(sigtab.D7.p.oi)

#Summarize sigtab.D7.p.oi
sum.sigtab.D7.p.oi <- summary(sigtab.D7.p.oi)
sum.sigtab.D7.p.oi

#ggplot
deseq.D7.p.oi <- ggplot(sigtab.D7.p.oi, aes(x=reorder(rownames(sigtab.D7.p.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.oi), y=-5e-07, label = paste(Phylum, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phylum\n in INFfeed Group Relative to INFnm\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', INFnm='#CC0066'))
deseq.D7.p.oi

#Add OTU and comparisons columns
sigtab.D7.p.oi
sigtab.D7.p.oi$OTU <- rownames(sigtab.D7.p.oi)
sigtab.D7.p.oi
sigtab.D7.p.oi$comp <- 'D7_INFfeed_vs_INFnm'

#Create final significant comparisons table
final.sigtab.phylum <- rbind(sigtab.D7.p.oi, final.sigtab.phylum)


######### 3. Day 7 INFfeed vs INFinject ###################

#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFinject")
#19
sum(meta$Set == "D7_INFfeed")
#13

#Extract results from a DESeq analysis, organize table
sample_data(FS9.D7.p)$Set <- factor(sample_data(FS9.D7.p)$Set,
                                  levels =c("D7_INFinject", "D7_INFfeed",
                                            'D7_INFnm'))
FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.p.De)
#[1] "Intercept"       "Set_D7_INFfeed_vs_D7_INFinject" "Set_D7_INFnm_vs_D7_INFinject"  
res.D7.p.oj = lfcShrink(FS9.D7.p.De, coef = "Set_D7_INFfeed_vs_D7_INFinject", type = 'apeglm')
sigtab.D7.p.oj = res.D7.p.oi[which(res.D7.p.oi$padj < .05), ]
sigtab.D7.p.oj = cbind(as(sigtab.D7.p.oj, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p.oj), ], "matrix"))
format(sigtab.D7.p.oj$padj, scientific = TRUE)
sigtab.D7.p.oj$newp <- format(round(sigtab.D7.p.oj$padj, digits = 3), scientific = TRUE)
sigtab.D7.p.oj$Treatment <- ifelse(sigtab.D7.p.oj$log2FoldChange >=0, "INFfeed", "INFinject")
head(sigtab.D7.p.oj)

#Summarize sigtab.D7.p.oj
sum.sigtab.D7.p.oj <- summary(sigtab.D7.p.oj)
sum.sigtab.D7.p.oj

#ggplot
deseq.D7.p.oj <- ggplot(sigtab.D7.p.oj, aes(x=reorder(rownames(sigtab.D7.p.oj), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.oj), y=0, label = paste(Phylum, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phylum\n in INFfeed Group Relative to INFinject\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', INFinject='#E69F00'))
deseq.D7.p.oj

#Add OTU and comparisons columns
sigtab.D7.p.oj
sigtab.D7.p.oj$OTU <- rownames(sigtab.D7.p.oj)
sigtab.D7.p.oj
sigtab.D7.p.oj$comp <- 'D7_INFfeed_vs_INFinject'

#Create final significant comparisons table
final.sigtab.phylum <- rbind(sigtab.D7.p.oj, final.sigtab.phylum)






#######################################################################################################

#write csv
write.csv(final.sigtab.phylum, file= "FS9_FinalDiffAbund_Phylum_OutDoubletons_Q2_D4D7.csv")

#######################################################################################################