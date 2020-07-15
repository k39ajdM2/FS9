#####################################################################################################
#FS9 DESeq2
#Kathy Mou

#Purpose: This code uses DESeq2 package to identify nasal microbial genera that were differentially 
#abundant between the four groups on the four days sampled

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
#R version 3.6.0 (2019-04-26)

########################################################################################################

#Import files
otu <- import_mothur(mothur_shared_file = './data/stability.outsingletons.abund.opti_mcc.shared') #use unrarified data
taxo <- import_mothur(mothur_constaxonomy_file = './data/stability.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
meta <- read.table(file = './data/FS9_metadata.csv', sep = ',', header = TRUE)

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
#otu_table()   OTU Table:         [ 2838 taxa and 172 samples ]
#sample_data() Sample Data:       [ 172 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 2838 taxa by 6 taxonomic ranks ]

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
# INFinject vs INFnm 
# INFfeed vs INFnm 
# INFfeed vs INFinject
# INFnm vs NONINFnm 
# INFinject vs NONINFnm
# INFfeed vs NONINFnm

# Day 0 
# INFinject vs INFnm 
# INFfeed vs INFnm 
# INFfeed vs INFinject
# INFnm vs NONINFnm 
# INFinject vs NONINFnm
# INFfeed vs NONINFnm

# Day 4 
# INFinject vs INFnm 
# INFfeed vs INFnm 
# INFfeed vs INFinject
# INFnm vs NONINFnm 
# INFinject vs NONINFnm
# INFfeed vs NONINFnm

# Day 7
# INFinject vs INFnm 
# INFfeed vs INFnm 
# INFfeed vs INFinject
# INFnm vs NONINFnm 
# INFinject vs NONINFnm
# INFfeed vs NONINFnm

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

######### 1. Day -3 INFinject vs INFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFnm")
#20
sum(meta$Set == "DNEG3_INFinject")
#20

resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                            "Set_DNEG3_INFinject_vs_DNEG3_INFfeed"
#[3] "Set_DNEG3_INFnm_vs_DNEG3_INFfeed"     "Set_DNEG3_NONINFnm_vs_DNEG3_INFfeed" 

#re-level your factor and re-run DESeq2
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c('DNEG3_INFnm',"DNEG3_NONINFnm",
                                               "DNEG3_INFinject", "DNEG3_INFfeed"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                          "Set_DNEG3_NONINFnm_vs_DNEG3_INFnm" 
#[3] "Set_DNEG3_INFinject_vs_DNEG3_INFnm" "Set_DNEG3_INFfeed_vs_DNEG3_INFnm"  
res.DNEG3.ji = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFinject_vs_DNEG3_INFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.DNEG3.ji = res.DNEG3.ji[which(res.DNEG3.ji$padj < .05), ]
sigtab.DNEG3.ji = cbind(as(sigtab.DNEG3.ji, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.ji), ], "matrix"))
sigtab.DNEG3.ji$newp <- format(round(sigtab.DNEG3.ji$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.ji$Treatment <- ifelse(sigtab.DNEG3.ji$log2FoldChange >=0, "INFinject", "INFnm")

deseq.DNEG3.ji <- 
  ggplot(sigtab.DNEG3.ji, aes(x=reorder(rownames(sigtab.DNEG3.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.ji), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFinject Group Relative to INFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066'))
deseq.DNEG3.ji

#Add OTU and comparisons columns
sigtab.DNEG3.ji
sigtab.DNEG3.ji$OTU <- rownames(sigtab.DNEG3.ji)
sigtab.DNEG3.ji
sigtab.DNEG3.ji$comp <- 'DNEG3_INFinject_vs_INFnm'

#Create final significant comparisons table
final.sigtab <- sigtab.DNEG3.ji



######### 2. Day -3 INFfeed vs INFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFnm")
#20
sum(meta$Set == "DNEG3_INFfeed")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                          "Set_DNEG3_NONINFnm_vs_DNEG3_INFnm" 
#[3] "Set_DNEG3_INFinject_vs_DNEG3_INFnm" "Set_DNEG3_INFfeed_vs_DNEG3_INFnm" 
res.DNEG3.oi = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFfeed_vs_DNEG3_INFnm", type = 'apeglm')
sigtab.DNEG3.oi = res.DNEG3.oi[which(res.DNEG3.oi$padj < .05), ]
sigtab.DNEG3.oi = cbind(as(sigtab.DNEG3.oi, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.oi), ], "matrix"))
format(sigtab.DNEG3.oi$padj, scientific = TRUE)
sigtab.DNEG3.oi$newp <- format(round(sigtab.DNEG3.oi$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.oi$Treatment <- ifelse(sigtab.DNEG3.oi$log2FoldChange >=0, "INFfeed", "INFnm")

#Summarize sigtab.DNEG3.oi
sum.sigtab.DNEG3.oi <- summary(sigtab.DNEG3.oi)
sum.sigtab.DNEG3.oi

#ggplot
deseq.DNEG3.oi <- ggplot(sigtab.DNEG3.oi, aes(x=reorder(rownames(sigtab.DNEG3.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.oi), y=-1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFfeed Group Relative to INFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', INFnm='#CC0066'))
deseq.DNEG3.oi

#Add OTU and comparisons columns
sigtab.DNEG3.oi
sigtab.DNEG3.oi$OTU <- rownames(sigtab.DNEG3.oi)
sigtab.DNEG3.oi
sigtab.DNEG3.oi$comp <- 'DNEG3_INFfeed_vs_INFnm'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.DNEG3.oi, final.sigtab)


######### 3. Day -3 INFfeed vs INFinject ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFinject")
#20
sum(meta$Set == "DNEG3_INFfeed")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"    
#[3] "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm" "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm"
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c("DNEG3_INFinject", "DNEG3_INFfeed",
                                               'DNEG3_INFnm',"DNEG3_NONINFnm"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFfeed_vs_DNEG3_INFinject" 
#[3] "Set_DNEG3_INFnm_vs_DNEG3_INFinject"    "Set_DNEG3_NONINFnm_vs_DNEG3_INFinject"
res.DNEG3.oj = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFfeed_vs_DNEG3_INFinject", type = 'apeglm')
sigtab.DNEG3.oj = res.DNEG3.oi[which(res.DNEG3.oi$padj < .05), ]
sigtab.DNEG3.oj = cbind(as(sigtab.DNEG3.oj, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.oj), ], "matrix"))
format(sigtab.DNEG3.oj$padj, scientific = TRUE)
sigtab.DNEG3.oj$newp <- format(round(sigtab.DNEG3.oj$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.oj$Treatment <- ifelse(sigtab.DNEG3.oj$log2FoldChange >=0, "INFfeed", "INFinject")

#Summarize sigtab.DNEG3.oj
sum.sigtab.DNEG3.oj <- summary(sigtab.DNEG3.oj)
sum.sigtab.DNEG3.oj

#ggplot
deseq.DNEG3.oj <- ggplot(sigtab.DNEG3.oj, aes(x=reorder(rownames(sigtab.DNEG3.oj), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.oj), y=-1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFfeed Group Relative to INFinject in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', INFinject='#E69F00'))
deseq.DNEG3.oj

#Add OTU and comparisons columns
sigtab.DNEG3.oj
sigtab.DNEG3.oj$OTU <- rownames(sigtab.DNEG3.oj)
sigtab.DNEG3.oj
sigtab.DNEG3.oj$comp <- 'DNEG3_INFfeed_vs_INFinject'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.DNEG3.oj, final.sigtab)





######### 4. Day -3 INFnm vs NONINFnm  ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_NONINFnm")
#20
sum(meta$Set == "DNEG3_INFnm")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c('DNEG3_NONINFnm','DNEG3_INFnm', 
                                               "DNEG3_INFinject", "DNEG3_INFfeed"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"    
#[3] "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm" "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm" 
res.DNEG3.in = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm", type = 'apeglm')
sigtab.DNEG3.in = res.DNEG3.in[which(res.DNEG3.in$padj < .05), ]
sigtab.DNEG3.in = cbind(as(sigtab.DNEG3.in, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.in), ], "matrix"))
format(sigtab.DNEG3.in$padj, scientific = TRUE)
sigtab.DNEG3.in$newp <- format(round(sigtab.DNEG3.in$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.in$Treatment <- ifelse(sigtab.DNEG3.in$log2FoldChange >=0, "INFnm", "NONINFnm")

#Summarize sigtab.DNEG3.in
sum.sigtab.DNEG3.in <- summary(sigtab.DNEG3.in)
sum.sigtab.DNEG3.in

#ggplot
deseq.DNEG3.in <- ggplot(sigtab.DNEG3.in, aes(x=reorder(rownames(sigtab.DNEG3.in), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.in), y=2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface= "italic")+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFnm Group Relative to NONINFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.DNEG3.in

#Add OTU and comparisons columns
sigtab.DNEG3.in
sigtab.DNEG3.in$OTU <- rownames(sigtab.DNEG3.in)
sigtab.DNEG3.in
sigtab.DNEG3.in$comp <- 'DNEG3_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.DNEG3.in)


######### 5. Day -3 INFinject vs NONINFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_NONINFnm")
#20
sum(meta$Set == "DNEG3_INFinject")
#20

resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFfeed_vs_DNEG3_INFinject" 
#[3] "Set_DNEG3_INFnm_vs_DNEG3_INFinject"    "Set_DNEG3_NONINFnm_vs_DNEG3_INFinject" 
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c("DNEG3_NONINFnm", 'DNEG3_INFnm',
                                               "DNEG3_INFinject", "DNEG3_INFfeed"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"    
#[3] "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm" "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm"  
res.DNEG3.jn = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm", type = 'apeglm')
sigtab.DNEG3.jn = res.DNEG3.jn[which(res.DNEG3.jn$padj < .05), ]
sigtab.DNEG3.jn = cbind(as(sigtab.DNEG3.jn, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.jn), ], "matrix"))
sigtab.DNEG3.jn$newp <- format(round(sigtab.DNEG3.jn$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.jn$Treatment <- ifelse(sigtab.DNEG3.jn$log2FoldChange >=0, "INFinject", "NONINFnm")

deseq.DNEG3.jn <- 
  ggplot(sigtab.DNEG3.jn, aes(x=reorder(rownames(sigtab.DNEG3.jn), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.jn), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFinject Group Relative to NONINFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', NONINFnm='#56B4E9'))
deseq.DNEG3.jn

#Add OTU and comparisons columns
sigtab.DNEG3.jn
sigtab.DNEG3.jn$OTU <- rownames(sigtab.DNEG3.jn)
sigtab.DNEG3.jn
sigtab.DNEG3.jn$comp <- 'DNEG3_INFinject_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.DNEG3.jn)


######### 6. Day -3 INFfeed vs NONINFnm  ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_NONINFnm")
#20
sum(meta$Set == "DNEG3_INFfeed")
#20

#Extract results from a DESeq analysis, organize table
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"    
#[3] "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm" "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm" 
res.DNEG3.on = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm", type = 'apeglm')
sigtab.DNEG3.on = res.DNEG3.on[which(res.DNEG3.on$padj < .05), ]
sigtab.DNEG3.on = cbind(as(sigtab.DNEG3.on, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.on), ], "matrix"))
format(sigtab.DNEG3.on$padj, scientific = TRUE)
sigtab.DNEG3.on$newp <- format(round(sigtab.DNEG3.on$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.on$Treatment <- ifelse(sigtab.DNEG3.on$log2FoldChange >=0, "INFfeed", "NONINFnm")

#Summarize sigtab.DNEG3.on
sum.sigtab.DNEG3.on <- summary(sigtab.DNEG3.on)
sum.sigtab.DNEG3.on

#ggplot
deseq.DNEG3.on <- ggplot(sigtab.DNEG3.on, aes(x=reorder(rownames(sigtab.DNEG3.on), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.on), y=-1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFfeed Group Relative to NONINFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', NONINFnm='#56B4E9'))
deseq.DNEG3.on

#Add OTU and comparisons columns
sigtab.DNEG3.on
sigtab.DNEG3.on$OTU <- rownames(sigtab.DNEG3.on)
sigtab.DNEG3.on
sigtab.DNEG3.on$comp <- 'DNEG3_INFfeed_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.DNEG3.on)


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

######### 1. Day 0 INFinject vs INFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_INFnm")
#8
sum(meta$Set == "D0_INFinject")
#13

resultsNames(FS9.D0.De)
#[1] [1] "Intercept"                      "Set_D0_INFinject_vs_D0_INFfeed" "Set_D0_INFnm_vs_D0_INFfeed"    
#[4] "Set_D0_NONINFnm_vs_D0_INFfeed"
sample_data(FS9.D0)$Set <- factor(sample_data(FS9.D0)$Set,
                                     levels =c('D0_INFnm',"D0_NONINFnm",
                                               "D0_INFinject", "D0_INFfeed"))
FS9.D0.De <- phyloseq_to_deseq2(FS9.D0, ~ Set)
FS9.D0.De <- DESeq(FS9.D0.De, test = "Wald", fitType = "parametric")
FS9.D0.De$Set
resultsNames(FS9.D0.De)
#[1] "Intercept"                    "Set_D0_NONINFnm_vs_D0_INFnm"  
#[3] "Set_D0_INFinject_vs_D0_INFnm" "Set_D0_INFfeed_vs_D0_INFnm"
res.D0.ji = lfcShrink(FS9.D0.De, coef = "Set_D0_INFinject_vs_D0_INFnm", type = 'apeglm')
sigtab.D0.ji = res.D0.ji[which(res.D0.ji$padj < .05), ]
sigtab.D0.ji = cbind(as(sigtab.D0.ji, "data.frame"), as(tax_table(FS9.D0)[rownames(sigtab.D0.ji), ], "matrix"))
sigtab.D0.ji$newp <- format(round(sigtab.D0.ji$padj, digits = 3), scientific = TRUE)
sigtab.D0.ji$Treatment <- ifelse(sigtab.D0.ji$log2FoldChange >=0, "INFinject", "INFnm")
head(sigtab.D0.ji) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison

#CONTINUE HERE
######### 2. Day 0 INFfeed vs INFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_INFnm")
#8
sum(meta$Set == "D0_INFfeed")
#5

#Extract results from a DESeq analysis, organize table
FS9.D0.De <- phyloseq_to_deseq2(FS9.D0, ~ Set)
FS9.D0.De <- DESeq(FS9.D0.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D0.De)
#1] "Intercept"                    "Set_D0_NONINFnm_vs_D0_INFnm"  "Set_D0_INFinject_vs_D0_INFnm"
#[4] "Set_D0_INFfeed_vs_D0_INFnm"   
res.D0.oi = lfcShrink(FS9.D0.De, coef = "Set_D0_INFfeed_vs_D0_INFnm", type = 'apeglm')
sigtab.D0.oi = res.D0.oi[which(res.D0.oi$padj < .05), ]
sigtab.D0.oi = cbind(as(sigtab.D0.oi, "data.frame"), as(tax_table(FS9.D0)[rownames(sigtab.D0.oi), ], "matrix"))
format(sigtab.D0.oi$padj, scientific = TRUE)
sigtab.D0.oi$newp <- format(round(sigtab.D0.oi$padj, digits = 3), scientific = TRUE)
sigtab.D0.oi$Treatment <- ifelse(sigtab.D0.oi$log2FoldChange >=0, "INFfeed", "INFnm")
head(sigtab.D0.oi) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


######### 3. Day 0 INFfeed vs INFinject ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_INFinject")
#13
sum(meta$Set == "D0_INFfeed")
#5

#Extract results from a DESeq analysis, organize table
FS9.D0.De$Set
resultsNames(FS9.D0.De)
#[1] "Intercept"                    "Set_D0_NONINFnm_vs_D0_INFnm"  "Set_D0_INFinject_vs_D0_INFnm"
#[4] "Set_D0_INFfeed_vs_D0_INFnm"  
sample_data(FS9.D0)$Set <- factor(sample_data(FS9.D0)$Set,
                                     levels =c("D0_INFinject", "D0_INFfeed",
                                               'D0_INFnm',"D0_NONINFnm"))
FS9.D0.De <- phyloseq_to_deseq2(FS9.D0, ~ Set)
FS9.D0.De <- DESeq(FS9.D0.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D0.De)
#[1] "Intercept"                       "Set_D0_INFfeed_vs_D0_INFinject"  "Set_D0_INFnm_vs_D0_INFinject"   
#[4] "Set_D0_NONINFnm_vs_D0_INFinject"
res.D0.oj = lfcShrink(FS9.D0.De, coef = "Set_D0_INFfeed_vs_D0_INFinject", type = 'apeglm')
sigtab.D0.oj = res.D0.oi[which(res.D0.oi$padj < .05), ]
sigtab.D0.oj = cbind(as(sigtab.D0.oj, "data.frame"), as(tax_table(FS9.D0)[rownames(sigtab.D0.oj), ], "matrix"))
format(sigtab.D0.oj$padj, scientific = TRUE)
sigtab.D0.oj$newp <- format(round(sigtab.D0.oj$padj, digits = 3), scientific = TRUE)
sigtab.D0.oj$Treatment <- ifelse(sigtab.D0.oj$log2FoldChange >=0, "INFfeed", "INFinject")
head(sigtab.D0.oj) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


######### 4. Day 0 INFnm vs NONINFnm  ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_NONINFnm")
#9
sum(meta$Set == "D0_INFnm")
#8

#Extract results from a DESeq analysis, organize table
resultsNames(FS9.D0.De)
#[1] "Intercept"                       "Set_D0_INFfeed_vs_D0_INFinject"  "Set_D0_INFnm_vs_D0_INFinject"   
#[4] "Set_D0_NONINFnm_vs_D0_INFinject"
sample_data(FS9.D0)$Set <- factor(sample_data(FS9.D0)$Set,
                                     levels =c('D0_NONINFnm','D0_INFnm', 
                                               "D0_INFinject", "D0_INFfeed"))
FS9.D0.De <- phyloseq_to_deseq2(FS9.D0, ~ Set)
FS9.D0.De <- DESeq(FS9.D0.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D0.De)
#[1] "Intercept"                             "Set_D0_INFnm_vs_D0_NONINFnm"    
#[3] "Set_D0_INFinject_vs_D0_NONINFnm" "Set_D0_INFfeed_vs_D0_NONINFnm" 
res.D0.in = lfcShrink(FS9.D0.De, coef = "Set_D0_INFnm_vs_D0_NONINFnm", type = 'apeglm')
sigtab.D0.in = res.D0.in[which(res.D0.in$padj < .05), ]
sigtab.D0.in = cbind(as(sigtab.D0.in, "data.frame"), as(tax_table(FS9.D0)[rownames(sigtab.D0.in), ], "matrix"))
format(sigtab.D0.in$padj, scientific = TRUE)
sigtab.D0.in$newp <- format(round(sigtab.D0.in$padj, digits = 3), scientific = TRUE)
sigtab.D0.in$Treatment <- ifelse(sigtab.D0.in$log2FoldChange >=0, "INFnm", "NONINFnm")

#Summarize sigtab.D0.in
sum.sigtab.D0.in <- summary(sigtab.D0.in)
sum.sigtab.D0.in

#ggplot
deseq.D0.in <- ggplot(sigtab.D0.in, aes(x=reorder(rownames(sigtab.D0.in), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.in), y=2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface= "italic")+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFnm Group Relative to NONINFnm in Fecal Microbiota on Day 0')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.D0.in

#Add OTU and comparisons columns
sigtab.D0.in
sigtab.D0.in$OTU <- rownames(sigtab.D0.in)
sigtab.D0.in
sigtab.D0.in$comp <- 'D0_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D0.in)


######### 5. Day 0 INFinject vs NONINFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_NONINFnm")
#9
sum(meta$Set == "D0_INFinject")
#13

resultsNames(FS9.D0.De)
#[1] "Intercept"                       "Set_D0_INFnm_vs_D0_NONINFnm"     "Set_D0_INFinject_vs_D0_NONINFnm"
#[4] "Set_D0_INFfeed_vs_D0_NONINFnm"   
res.D0.jn = lfcShrink(FS9.D0.De, coef = "Set_D0_INFinject_vs_D0_NONINFnm", type = 'apeglm')
sigtab.D0.jn = res.D0.jn[which(res.D0.jn$padj < .05), ]
sigtab.D0.jn = cbind(as(sigtab.D0.jn, "data.frame"), as(tax_table(FS9.D0)[rownames(sigtab.D0.jn), ], "matrix"))
sigtab.D0.jn$newp <- format(round(sigtab.D0.jn$padj, digits = 3), scientific = TRUE)
sigtab.D0.jn$Treatment <- ifelse(sigtab.D0.jn$log2FoldChange >=0, "INFinject", "NONINFnm")
head(sigtab.D0.jn)

deseq.D0.jn <- 
  ggplot(sigtab.D0.jn, aes(x=reorder(rownames(sigtab.D0.jn), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.jn), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFinject Group Relative to NONINFnm in Fecal Microbiota on Day 0')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', NONINFnm='#56B4E9'))
deseq.D0.jn

#Add OTU and comparisons columns
sigtab.D0.jn
sigtab.D0.jn$OTU <- rownames(sigtab.D0.jn)
sigtab.D0.jn
sigtab.D0.jn$comp <- 'D0_INFinject_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D0.jn)


######### 6. Day 0 INFfeed vs NONINFnm  ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_NONINFnm")
#9
sum(meta$Set == "D0_INFfeed")
#5

#Extract results from a DESeq analysis, organize table
resultsNames(FS9.D0.De)
#[1] "Intercept"                             "Set_D0_INFnm_vs_D0_NONINFnm"    
#[3] "Set_D0_INFinject_vs_D0_NONINFnm" "Set_D0_INFfeed_vs_D0_NONINFnm" 
res.D0.on = lfcShrink(FS9.D0.De, coef = "Set_D0_INFfeed_vs_D0_NONINFnm", type = 'apeglm')
sigtab.D0.on = res.D0.on[which(res.D0.on$padj < .05), ]
sigtab.D0.on = cbind(as(sigtab.D0.on, "data.frame"), as(tax_table(FS9.D0)[rownames(sigtab.D0.on), ], "matrix"))
format(sigtab.D0.on$padj, scientific = TRUE)
sigtab.D0.on$newp <- format(round(sigtab.D0.on$padj, digits = 3), scientific = TRUE)
sigtab.D0.on$Treatment <- ifelse(sigtab.D0.on$log2FoldChange >=0, "INFfeed", "NONINFnm")
head(sigtab.D0.on) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


##################################################### Day 4 ######################################################################

sample_data(FS9.order)

FS9.D4 <- subset_samples(FS9.order, Day == 'D4')
sample_sums(FS9.D4)
colnames(otu_table(FS9.D4)) #check on all the sample names
FS9.D4 <- prune_taxa(taxa_sums(FS9.D4) > 1, FS9.D4)
#if taxa_sums is >1, then it will print that out in FS9.D4 object and not include anything with <1.
rowSums(FS9.D4@otu_table)

#Look at what Set is
sample_data(FS9.D4)
FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")

######### 1. Day 4 INFinject vs INFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_INFnm")
#6
sum(meta$Set == "D4_INFinject")
#7

resultsNames(FS9.D4.De)
#[1] "Intercept"                      "Set_D4_INFinject_vs_D4_INFfeed" 
#[3] "Set_D4_INFnm_vs_D4_INFfeed"     "Set_D4_NONINFnm_vs_D4_INFfeed" 
sample_data(FS9.D4)$Set <- factor(sample_data(FS9.D4)$Set,
                                     levels =c('D4_INFnm',"D4_NONINFnm",
                                               "D4_INFinject", "D4_INFfeed"))
FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")
FS9.D4.De$Set
resultsNames(FS9.D4.De)
#[1] "Intercept"                    "Set_D4_NONINFnm_vs_D4_INFnm"  
#[3] "Set_D4_INFinject_vs_D4_INFnm" "Set_D4_INFfeed_vs_D4_INFnm" 
res.D4.ji = lfcShrink(FS9.D4.De, coef = "Set_D4_INFinject_vs_D4_INFnm", type = 'apeglm')
sigtab.D4.ji = res.D4.ji[which(res.D4.ji$padj < .05), ]
sigtab.D4.ji = cbind(as(sigtab.D4.ji, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.ji), ], "matrix"))
sigtab.D4.ji$newp <- format(round(sigtab.D4.ji$padj, digits = 3), scientific = TRUE)
sigtab.D4.ji$Treatment <- ifelse(sigtab.D4.ji$log2FoldChange >=0, "INFinject", "INFnm")
head(sigtab.D4.ji) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


######### 2. Day 4 INFfeed vs INFnm ###################

#NONINFnm = N
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
#[1] "Intercept"                    "Set_D4_NONINFnm_vs_D4_INFnm"  
#[3] "Set_D4_INFinject_vs_D4_INFnm" "Set_D4_INFfeed_vs_D4_INFnm" 
res.D4.oi = lfcShrink(FS9.D4.De, coef = "Set_D4_INFfeed_vs_D4_INFnm", type = 'apeglm')
sigtab.D4.oi = res.D4.oi[which(res.D4.oi$padj < .05), ]
sigtab.D4.oi = cbind(as(sigtab.D4.oi, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.oi), ], "matrix"))
format(sigtab.D4.oi$padj, scientific = TRUE)
sigtab.D4.oi$newp <- format(round(sigtab.D4.oi$padj, digits = 3), scientific = TRUE)
sigtab.D4.oi$Treatment <- ifelse(sigtab.D4.oi$log2FoldChange >=0, "INFfeed", "INFnm")
head(sigtab.D4.oi)  #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


######### 3. Day 4 INFfeed vs INFinject ###################

#NONINFnm = N
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
#[1] "Intercept"                    "Set_D4_NONINFnm_vs_D4_INFnm"  
#[3] "Set_D4_INFinject_vs_D4_INFnm" "Set_D4_INFfeed_vs_D4_INFnm"
sample_data(FS9.D4)$Set <- factor(sample_data(FS9.D4)$Set,
                                     levels =c("D4_INFinject", "D4_INFfeed",
                                               'D4_INFnm',"D4_NONINFnm"))
FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D4.De)
#[1] "Intercept"                       "Set_D4_INFfeed_vs_D4_INFinject"  "Set_D4_INFnm_vs_D4_INFinject"   
#[4] "Set_D4_NONINFnm_vs_D4_INFinject"
res.D4.oj = lfcShrink(FS9.D4.De, coef = "Set_D4_INFfeed_vs_D4_INFinject", type = 'apeglm')
sigtab.D4.oj = res.D4.oi[which(res.D4.oi$padj < .05), ]
sigtab.D4.oj = cbind(as(sigtab.D4.oj, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.oj), ], "matrix"))
format(sigtab.D4.oj$padj, scientific = TRUE)
sigtab.D4.oj$newp <- format(round(sigtab.D4.oj$padj, digits = 3), scientific = TRUE)
sigtab.D4.oj$Treatment <- ifelse(sigtab.D4.oj$log2FoldChange >=0, "INFfeed", "INFinject")
head(sigtab.D4.oj) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
#in abundance between the two groups, so I will skip to the next comparison


######### 4. Day 4 INFnm vs NONINFnm  ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_NONINFnm")
#5
sum(meta$Set == "D4_INFnm")
#6

#Extract results from a DESeq analysis, organize table
FS9.D4.De$Set
resultsNames(FS9.D4.De)
#[1] "Intercept"                       "Set_D4_INFfeed_vs_D4_INFinject"  "Set_D4_INFnm_vs_D4_INFinject"   
#[4] "Set_D4_NONINFnm_vs_D4_INFinject"
sample_data(FS9.D4)$Set <- factor(sample_data(FS9.D4)$Set,
                                     levels =c('D4_NONINFnm','D4_INFnm', 
                                               "D4_INFinject", "D4_INFfeed"))
FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D4.De)
#[1] "Intercept"                       "Set_D4_INFnm_vs_D4_NONINFnm"     "Set_D4_INFinject_vs_D4_NONINFnm"
#[4] "Set_D4_INFfeed_vs_D4_NONINFnm"
res.D4.in = lfcShrink(FS9.D4.De, coef = "Set_D4_INFnm_vs_D4_NONINFnm", type = 'apeglm')
sigtab.D4.in = res.D4.in[which(res.D4.in$padj < .05), ]
sigtab.D4.in = cbind(as(sigtab.D4.in, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.in), ], "matrix"))
format(sigtab.D4.in$padj, scientific = TRUE)
sigtab.D4.in$newp <- format(round(sigtab.D4.in$padj, digits = 3), scientific = TRUE)
sigtab.D4.in$Treatment <- ifelse(sigtab.D4.in$log2FoldChange >=0, "INFnm", "NONINFnm")

#Summarize sigtab.D4.in
sum.sigtab.D4.in <- summary(sigtab.D4.in)
sum.sigtab.D4.in

#ggplot
deseq.D4.in <- ggplot(sigtab.D4.in, aes(x=reorder(rownames(sigtab.D4.in), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.in), y=2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface= "italic")+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFnm Group Relative to NONINFnm in Fecal Microbiota on Day 4')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.D4.in

#Add OTU and comparisons columns
sigtab.D4.in
sigtab.D4.in$OTU <- rownames(sigtab.D4.in)
sigtab.D4.in
sigtab.D4.in$comp <- 'D4_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D4.in)


######### 5. Day 4 INFinject vs NONINFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_NONINFnm")
#5
sum(meta$Set == "D4_INFinject")
#7

resultsNames(FS9.D4.De)
#[1] "Intercept"                       "Set_D4_INFnm_vs_D4_NONINFnm"     "Set_D4_INFinject_vs_D4_NONINFnm"
#[4] "Set_D4_INFfeed_vs_D4_NONINFnm" 
sample_data(FS9.D4)$Set <- factor(sample_data(FS9.D4)$Set,
                                     levels =c("D4_NONINFnm", 'D4_INFnm',
                                               "D4_INFinject", "D4_INFfeed"))
FS9.D4.De <- phyloseq_to_deseq2(FS9.D4, ~ Set)
FS9.D4.De <- DESeq(FS9.D4.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D4.De)
#[1] "Intercept"                       "Set_D4_INFnm_vs_D4_NONINFnm"     "Set_D4_INFinject_vs_D4_NONINFnm"
#[4] "Set_D4_INFfeed_vs_D4_NONINFnm"  
res.D4.jn = lfcShrink(FS9.D4.De, coef = "Set_D4_INFinject_vs_D4_NONINFnm", type = 'apeglm')
sigtab.D4.jn = res.D4.jn[which(res.D4.jn$padj < .05), ]
sigtab.D4.jn = cbind(as(sigtab.D4.jn, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.jn), ], "matrix"))
sigtab.D4.jn$newp <- format(round(sigtab.D4.jn$padj, digits = 3), scientific = TRUE)
sigtab.D4.jn$Treatment <- ifelse(sigtab.D4.jn$log2FoldChange >=0, "INFinject", "NONINFnm")

deseq.D4.jn <- 
  ggplot(sigtab.D4.jn, aes(x=reorder(rownames(sigtab.D4.jn), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.jn), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFinject Group Relative to NONINFnm in Fecal Microbiota on Day 4')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', NONINFnm='#56B4E9'))
deseq.D4.jn

#Add OTU and comparisons columns
sigtab.D4.jn
sigtab.D4.jn$OTU <- rownames(sigtab.D4.jn)
sigtab.D4.jn
sigtab.D4.jn$comp <- 'D4_INFinject_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D4.jn)


######### 6. Day 4 INFfeed vs NONINFnm  ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_NONINFnm")
#5
sum(meta$Set == "D4_INFfeed")
#4

#Extract results from a DESeq analysis, organize table
resultsNames(FS9.D4.De)
#[1] "Intercept"                       "Set_D4_INFnm_vs_D4_NONINFnm"     "Set_D4_INFinject_vs_D4_NONINFnm"
#[4] "Set_D4_INFfeed_vs_D4_NONINFnm"  
res.D4.on = lfcShrink(FS9.D4.De, coef = "Set_D4_INFfeed_vs_D4_NONINFnm", type = 'apeglm')
sigtab.D4.on = res.D4.on[which(res.D4.on$padj < .05), ]
sigtab.D4.on = cbind(as(sigtab.D4.on, "data.frame"), as(tax_table(FS9.D4)[rownames(sigtab.D4.on), ], "matrix"))
format(sigtab.D4.on$padj, scientific = TRUE)
sigtab.D4.on$newp <- format(round(sigtab.D4.on$padj, digits = 3), scientific = TRUE)
sigtab.D4.on$Treatment <- ifelse(sigtab.D4.on$log2FoldChange >=0, "INFfeed", "NONINFnm")
head(sigtab.D4.on) #DataFrame with 0 rows and 7 columns, meaning there were no orders that were significantly different
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

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFnm")
#20
sum(meta$Set == "DNEG3_INFinject")
#20

resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                            "Set_DNEG3_INFinject_vs_DNEG3_INFfeed"
#[3] "Set_DNEG3_INFnm_vs_DNEG3_INFfeed"     "Set_DNEG3_NONINFnm_vs_DNEG3_INFfeed" 

#re-level your factor and re-run DESeq2
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c('DNEG3_INFnm',"DNEG3_NONINFnm",
                                               "DNEG3_INFinject", "DNEG3_INFfeed"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                          "Set_DNEG3_NONINFnm_vs_DNEG3_INFnm" 
#[3] "Set_DNEG3_INFinject_vs_DNEG3_INFnm" "Set_DNEG3_INFfeed_vs_DNEG3_INFnm"  
res.DNEG3.ji = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFinject_vs_DNEG3_INFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (.3_INF_InjOTC)
# negative log2foldchanges are associated with the second group from this line (.3_INF_NoTRMT)
sigtab.DNEG3.ji = res.DNEG3.ji[which(res.DNEG3.ji$padj < .05), ]
sigtab.DNEG3.ji = cbind(as(sigtab.DNEG3.ji, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.ji), ], "matrix"))
sigtab.DNEG3.ji$newp <- format(round(sigtab.DNEG3.ji$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.ji$Treatment <- ifelse(sigtab.DNEG3.ji$log2FoldChange >=0, "INFinject", "INFnm")

deseq.DNEG3.ji <- 
  ggplot(sigtab.DNEG3.ji, aes(x=reorder(rownames(sigtab.DNEG3.ji), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.ji), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFinject Group Relative to INFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', INFnm='#CC0066'))
deseq.DNEG3.ji

#Add OTU and comparisons columns
sigtab.DNEG3.ji
sigtab.DNEG3.ji$OTU <- rownames(sigtab.DNEG3.ji)
sigtab.DNEG3.ji
sigtab.DNEG3.ji$comp <- 'DNEG3_INFinject_vs_INFnm'

#Create final significant comparisons table
final.sigtab <- sigtab.DNEG3.ji



######### 2. Day -3 INFfeed vs INFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFnm")
#20
sum(meta$Set == "DNEG3_INFfeed")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                          "Set_DNEG3_NONINFnm_vs_DNEG3_INFnm" 
#[3] "Set_DNEG3_INFinject_vs_DNEG3_INFnm" "Set_DNEG3_INFfeed_vs_DNEG3_INFnm" 
res.DNEG3.oi = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFfeed_vs_DNEG3_INFnm", type = 'apeglm')
sigtab.DNEG3.oi = res.DNEG3.oi[which(res.DNEG3.oi$padj < .05), ]
sigtab.DNEG3.oi = cbind(as(sigtab.DNEG3.oi, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.oi), ], "matrix"))
format(sigtab.DNEG3.oi$padj, scientific = TRUE)
sigtab.DNEG3.oi$newp <- format(round(sigtab.DNEG3.oi$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.oi$Treatment <- ifelse(sigtab.DNEG3.oi$log2FoldChange >=0, "INFfeed", "INFnm")

#Summarize sigtab.DNEG3.oi
sum.sigtab.DNEG3.oi <- summary(sigtab.DNEG3.oi)
sum.sigtab.DNEG3.oi

#ggplot
deseq.DNEG3.oi <- ggplot(sigtab.DNEG3.oi, aes(x=reorder(rownames(sigtab.DNEG3.oi), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.oi), y=-1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFfeed Group Relative to INFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', INFnm='#CC0066'))
deseq.DNEG3.oi

#Add OTU and comparisons columns
sigtab.DNEG3.oi
sigtab.DNEG3.oi$OTU <- rownames(sigtab.DNEG3.oi)
sigtab.DNEG3.oi
sigtab.DNEG3.oi$comp <- 'DNEG3_INFfeed_vs_INFnm'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.DNEG3.oi, final.sigtab)


######### 3. Day -3 INFfeed vs INFinject ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFinject")
#20
sum(meta$Set == "DNEG3_INFfeed")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"    
#[3] "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm" "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm"
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c("DNEG3_INFinject", "DNEG3_INFfeed",
                                               'DNEG3_INFnm',"DNEG3_NONINFnm"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFfeed_vs_DNEG3_INFinject" 
#[3] "Set_DNEG3_INFnm_vs_DNEG3_INFinject"    "Set_DNEG3_NONINFnm_vs_DNEG3_INFinject"
res.DNEG3.oj = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFfeed_vs_DNEG3_INFinject", type = 'apeglm')
sigtab.DNEG3.oj = res.DNEG3.oi[which(res.DNEG3.oi$padj < .05), ]
sigtab.DNEG3.oj = cbind(as(sigtab.DNEG3.oj, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.oj), ], "matrix"))
format(sigtab.DNEG3.oj$padj, scientific = TRUE)
sigtab.DNEG3.oj$newp <- format(round(sigtab.DNEG3.oj$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.oj$Treatment <- ifelse(sigtab.DNEG3.oj$log2FoldChange >=0, "INFfeed", "INFinject")

#Summarize sigtab.DNEG3.oj
sum.sigtab.DNEG3.oj <- summary(sigtab.DNEG3.oj)
sum.sigtab.DNEG3.oj

#ggplot
deseq.DNEG3.oj <- ggplot(sigtab.DNEG3.oj, aes(x=reorder(rownames(sigtab.DNEG3.oj), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.oj), y=-1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFfeed Group Relative to INFinject in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', INFinject='#E69F00'))
deseq.DNEG3.oj

#Add OTU and comparisons columns
sigtab.DNEG3.oj
sigtab.DNEG3.oj$OTU <- rownames(sigtab.DNEG3.oj)
sigtab.DNEG3.oj
sigtab.DNEG3.oj$comp <- 'DNEG3_INFfeed_vs_INFinject'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.DNEG3.oj, final.sigtab)





######### 4. Day -3 INFnm vs NONINFnm  ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O


meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_NONINFnm")
#20
sum(meta$Set == "DNEG3_INFnm")
#20

#Extract results from a DESeq analysis, organize table
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c('DNEG3_NONINFnm','DNEG3_INFnm', 
                                               "DNEG3_INFinject", "DNEG3_INFfeed"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"    
#[3] "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm" "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm" 
res.DNEG3.in = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm", type = 'apeglm')
sigtab.DNEG3.in = res.DNEG3.in[which(res.DNEG3.in$padj < .05), ]
sigtab.DNEG3.in = cbind(as(sigtab.DNEG3.in, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.in), ], "matrix"))
format(sigtab.DNEG3.in$padj, scientific = TRUE)
sigtab.DNEG3.in$newp <- format(round(sigtab.DNEG3.in$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.in$Treatment <- ifelse(sigtab.DNEG3.in$log2FoldChange >=0, "INFnm", "NONINFnm")

#Summarize sigtab.DNEG3.in
sum.sigtab.DNEG3.in <- summary(sigtab.DNEG3.in)
sum.sigtab.DNEG3.in

#ggplot
deseq.DNEG3.in <- ggplot(sigtab.DNEG3.in, aes(x=reorder(rownames(sigtab.DNEG3.in), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.in), y=2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface= "italic")+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFnm Group Relative to NONINFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.DNEG3.in

#Add OTU and comparisons columns
sigtab.DNEG3.in
sigtab.DNEG3.in$OTU <- rownames(sigtab.DNEG3.in)
sigtab.DNEG3.in
sigtab.DNEG3.in$comp <- 'DNEG3_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.DNEG3.in)


######### 5. Day -3 INFinject vs NONINFnm ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_NONINFnm")
#20
sum(meta$Set == "DNEG3_INFinject")
#20

resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFfeed_vs_DNEG3_INFinject" 
#[3] "Set_DNEG3_INFnm_vs_DNEG3_INFinject"    "Set_DNEG3_NONINFnm_vs_DNEG3_INFinject" 
sample_data(FS9.DNEG3)$Set <- factor(sample_data(FS9.DNEG3)$Set,
                                     levels =c("DNEG3_NONINFnm", 'DNEG3_INFnm',
                                               "DNEG3_INFinject", "DNEG3_INFfeed"))
FS9.DNEG3.De <- phyloseq_to_deseq2(FS9.DNEG3, ~ Set)
FS9.DNEG3.De <- DESeq(FS9.DNEG3.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.De$Set
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"    
#[3] "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm" "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm"  
res.DNEG3.jn = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm", type = 'apeglm')
sigtab.DNEG3.jn = res.DNEG3.jn[which(res.DNEG3.jn$padj < .05), ]
sigtab.DNEG3.jn = cbind(as(sigtab.DNEG3.jn, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.jn), ], "matrix"))
sigtab.DNEG3.jn$newp <- format(round(sigtab.DNEG3.jn$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.jn$Treatment <- ifelse(sigtab.DNEG3.jn$log2FoldChange >=0, "INFinject", "NONINFnm")

deseq.DNEG3.jn <- 
  ggplot(sigtab.DNEG3.jn, aes(x=reorder(rownames(sigtab.DNEG3.jn), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.jn), y=-2, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFinject Group Relative to NONINFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFinject='#E69F00', NONINFnm='#56B4E9'))
deseq.DNEG3.jn

#Add OTU and comparisons columns
sigtab.DNEG3.jn
sigtab.DNEG3.jn$OTU <- rownames(sigtab.DNEG3.jn)
sigtab.DNEG3.jn
sigtab.DNEG3.jn$comp <- 'DNEG3_INFinject_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.DNEG3.jn)


######### 6. Day -3 INFfeed vs NONINFnm  ###################

#NONINFnm = N
#INFnm = I
#INFinject= J
#INFfeed = O

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_NONINFnm")
#20
sum(meta$Set == "DNEG3_INFfeed")
#20

#Extract results from a DESeq analysis, organize table
resultsNames(FS9.DNEG3.De)
#[1] "Intercept"                             "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"    
#[3] "Set_DNEG3_INFinject_vs_DNEG3_NONINFnm" "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm" 
res.DNEG3.on = lfcShrink(FS9.DNEG3.De, coef = "Set_DNEG3_INFfeed_vs_DNEG3_NONINFnm", type = 'apeglm')
sigtab.DNEG3.on = res.DNEG3.on[which(res.DNEG3.on$padj < .05), ]
sigtab.DNEG3.on = cbind(as(sigtab.DNEG3.on, "data.frame"), as(tax_table(FS9.DNEG3)[rownames(sigtab.DNEG3.on), ], "matrix"))
format(sigtab.DNEG3.on$padj, scientific = TRUE)
sigtab.DNEG3.on$newp <- format(round(sigtab.DNEG3.on$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.on$Treatment <- ifelse(sigtab.DNEG3.on$log2FoldChange >=0, "INFfeed", "NONINFnm")

#Summarize sigtab.DNEG3.on
sum.sigtab.DNEG3.on <- summary(sigtab.DNEG3.on)
sum.sigtab.DNEG3.on

#ggplot
deseq.DNEG3.on <- ggplot(sigtab.DNEG3.on, aes(x=reorder(rownames(sigtab.DNEG3.on), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.on), y=-1, label = paste(Phylum,Order, sep = ' ')), size=5, fontface = 'italic')+ labs(x="Phylum Order")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant OTUs in INFfeed Group Relative to NONINFnm in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFfeed='#999999', NONINFnm='#56B4E9'))
deseq.DNEG3.on

#Add OTU and comparisons columns
sigtab.DNEG3.on
sigtab.DNEG3.on$OTU <- rownames(sigtab.DNEG3.on)
sigtab.DNEG3.on
sigtab.DNEG3.on$comp <- 'DNEG3_INFfeed_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.DNEG3.on)



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
sample_data(FS9.DNEG3.p)$Set <- factor(sample_data(FS9.DNEG3.p)$Set,
                                     levels =c('-3_INF_NoTRMT',"-3_NONINF_NoTRMT",
                                               "-3_INF_InjOTC", "-3_INF_OralOTC"))
FS9.DNEG3.p.De <- phyloseq_to_deseq2(FS9.DNEG3.p, ~ Set)
FS9.DNEG3.p.De <- DESeq(FS9.DNEG3.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.DNEG3.p.De)
res.DNEG3.p.ji = lfcShrink(FS9.DNEG3.p.De, coef = "Set_.3_INF_InjOTC_vs_.3_INF_NoTRMT", type = 'apeglm')
#res.DNEG3.p.ji = results(FS9.DNEG3.p.De, contrast=c("Set", "-3_INF_InjOTC", "-3_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
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
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p.ji), y=0, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
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
sample_data(FS9.DNEG3.p)$Set <- factor(sample_data(FS9.DNEG3.p)$Set,
                                       levels =c('-3_INF_NoTRMT',"-3_NONINF_NoTRMT",
                                                 "-3_INF_InjOTC", "-3_INF_OralOTC"))
FS9.DNEG3.p.De <- phyloseq_to_deseq2(FS9.DNEG3.p, ~ Set)
FS9.DNEG3.p.De <- DESeq(FS9.DNEG3.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.DNEG3.p.De)
#[1] "Intercept"                             "Set_.3_NONINF_NoTRMT_vs_.3_INF_NoTRMT" "Set_.3_INF_InjOTC_vs_.3_INF_NoTRMT"   
#[4] "Set_.3_INF_OralOTC_vs_.3_INF_NoTRMT"  
res.DNEG3.p.oi = lfcShrink(FS9.DNEG3.p.De, coef = "Set_.3_INF_OralOTC_vs_.3_INF_NoTRMT", type = 'apeglm')
#res.DNEG3.p.oi = results(FS9.DNEG3.p.De, contrast=c("Set", "-3_INF_OralOTC", "-3_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS9.DNEG3.p.De, contrast=c("Set","-3_INF_OralOTC", "-3_INF_NoTRMT")) 
sigtab.DNEG3.p.oi = res.DNEG3.p.oi[which(res.DNEG3.p.oi$padj < .05), ]
sigtab.DNEG3.p.oi = cbind(as(sigtab.DNEG3.p.oi, "data.frame"), as(tax_table(FS9.DNEG3.p)[rownames(sigtab.DNEG3.p.oi), ], "matrix"))
#Error in dimnames(x) <- dn : 
#length of 'dimnames' [1] not equal to array extent

sigtab.DNEG3.p.oi
#log2 fold change (MLE): Set -3_INF_OralOTC vs -3_INF_NoTRMT 
#Wald test p-value: Set -3_INF_OralOTC vs -3_INF_NoTRMT 
#DataFrame with 0 rows and 5 columns
#0 rows and 5 columns indicate there were no phyla that were found significantly differentially abundant between INF_OralOTC and INF_NoTRMT.
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
sample_data(FS9.DNEG3.p)$Set <- factor(sample_data(FS9.DNEG3.p)$Set,
                                       levels =c("-3_NONINF_NoTRMT",'-3_INF_NoTRMT',
                                                 "-3_INF_InjOTC", "-3_INF_OralOTC"))
FS9.DNEG3.p.De <- phyloseq_to_deseq2(FS9.DNEG3.p, ~ Set)
FS9.DNEG3.p.De <- DESeq(FS9.DNEG3.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.DNEG3.p.De)
#[1] "Intercept"                              "Set_.3_INF_NoTRMT_vs_.3_NONINF_NoTRMT"  "Set_.3_INF_InjOTC_vs_.3_NONINF_NoTRMT" 
#[4] "Set_.3_INF_OralOTC_vs_.3_NONINF_NoTRMT"  
res.DNEG3.p.in = lfcShrink(FS9.DNEG3.p.De, coef = "Set_.3_INF_NoTRMT_vs_.3_NONINF_NoTRMT", type = 'apeglm')
#res.DNEG3.p.in = results(FS9.DNEG3.p.De, contrast=c("Set", "-3_INF_NoTRMT", "-3_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
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
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p.in), y=0, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
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
sample_data(FS9.DNEG3.p)$Set <- factor(sample_data(FS9.DNEG3.p)$Set,
                                       levels =c("-3_NONINF_NoTRMT",'-3_INF_NoTRMT',
                                                 "-3_INF_InjOTC", "-3_INF_OralOTC"))
FS9.DNEG3.p.De <- phyloseq_to_deseq2(FS9.DNEG3.p, ~ Set)
FS9.DNEG3.p.De <- DESeq(FS9.DNEG3.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.DNEG3.p.De)
#[1] "Intercept"                              "Set_.3_INF_NoTRMT_vs_.3_NONINF_NoTRMT"  "Set_.3_INF_InjOTC_vs_.3_NONINF_NoTRMT" 
#[4] "Set_.3_INF_OralOTC_vs_.3_NONINF_NoTRMT"  
res.DNEG3.p.on = lfcShrink(FS9.DNEG3.p.De, coef = "Set_.3_INF_OralOTC_vs_.3_NONINF_NoTRMT", type = 'apeglm')
#res.DNEG3.p.on = results(FS9.DNEG3.p.De, contrast=c("Set", "-3_INF_OralOTC", "-3_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
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
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p.on), y=0, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
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


##################################################### Day 4 ######################################################################

sample_data(FS9.phylum)

FS9.D4.p <- subset_samples(FS9.phylum, Day == '4')
sample_sums(FS9.D4.p)
colnames(otu_table(FS9.D4.p)) #check on all the sample names
FS9.D4.p <- prune_taxa(taxa_sums(FS9.D4.p) > 1, FS9.D4.p)
#if taxa_sums is >1, then it will print that out in FS9.D4.p object and not include anything with <1.
rowSums(FS9.D4.p@otu_table)

#Look at what Set is
sample_data(FS9.D4.p)

FS9.D4.p.De <- phyloseq_to_deseq2(FS9.D4.p, ~ Set)
FS9.D4.p.De <- DESeq(FS9.D4.p.De, test = "Wald", fitType = "parametric")





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
sample_data(FS9.D7.p)$Set <- factor(sample_data(FS9.D7.p)$Set,
                                       levels =c('7_INF_NoTRMT', "7_NONINF_NoTRMT",
                                                 "7_INF_InjOTC", "7_INF_OralOTC"))
FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.p.De)
#[1] "Intercept"                           "Set_7_NONINF_NoTRMT_vs_7_INF_NoTRMT" "Set_7_INF_InjOTC_vs_7_INF_NoTRMT"   
#[4] "Set_7_INF_OralOTC_vs_7_INF_NoTRMT" 
res.D7.p.ji = lfcShrink(FS9.D7.p.De, coef = "Set_7_INF_InjOTC_vs_7_INF_NoTRMT", type = 'apeglm')
#res.D7.p.ji = results(FS9.D7.p.De, contrast=c("Set", "7_INF_InjOTC", "7_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
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
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.ji), y=0, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
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
sample_data(FS9.D7.p)$Set <- factor(sample_data(FS9.D7.p)$Set,
                                    levels =c('7_INF_NoTRMT', "7_NONINF_NoTRMT",
                                              "7_INF_InjOTC", "7_INF_OralOTC"))
FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.p.De)
#[1] "Intercept"                           "Set_7_NONINF_NoTRMT_vs_7_INF_NoTRMT" "Set_7_INF_InjOTC_vs_7_INF_NoTRMT"   
#[4] "Set_7_INF_OralOTC_vs_7_INF_NoTRMT"  
res.D7.p.oi = lfcShrink(FS9.D7.p.De, coef = "Set_7_INF_OralOTC_vs_7_INF_NoTRMT", type = 'apeglm')
#res.D7.p.oi = results(FS9.D7.p.De, contrast=c("Set", "7_INF_OralOTC", "7_INF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
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
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.oi), y=0, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
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
sample_data(FS9.D7.p)$Set <- factor(sample_data(FS9.D7.p)$Set,
                                    levels =c("7_INF_InjOTC", "7_INF_OralOTC",
                                              '7_INF_NoTRMT', "7_NONINF_NoTRMT"))
FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.p.De)
#[1] "Intercept"                           "Set_7_INF_OralOTC_vs_7_INF_InjOTC"   "Set_7_INF_NoTRMT_vs_7_INF_InjOTC"   
#[4] "Set_7_NONINF_NoTRMT_vs_7_INF_InjOTC"
res.D7.p.oj = lfcShrink(FS9.D7.p.De, coef = "Set_7_INF_OralOTC_vs_7_INF_InjOTC", type = 'apeglm')
#res.D7.p.oj = results(FS9.D7.p.De, contrast=c("Set", "7_INF_OralOTC", "7_INF_InjOTC"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
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
sample_data(FS9.D7.p)$Set <- factor(sample_data(FS9.D7.p)$Set,
                                    levels =c("7_NONINF_NoTRMT", '7_INF_NoTRMT', 
                                              "7_INF_InjOTC", "7_INF_OralOTC"))
FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.p.De)
#[1] "Intercept"                            "Set_7_INF_NoTRMT_vs_7_NONINF_NoTRMT"  "Set_7_INF_InjOTC_vs_7_NONINF_NoTRMT" 
#[4] "Set_7_INF_OralOTC_vs_7_NONINF_NoTRMT"  
res.D7.p.jn = lfcShrink(FS9.D7.p.De, coef = "Set_7_INF_InjOTC_vs_7_NONINF_NoTRMT", type = 'apeglm')
#res.D7.p.jn = results(FS9.D7.p.De, contrast=c("Set", "7_INF_InjOTC", "7_NONINF_NoTRMT"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
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
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p.jn), y=0, label = paste(Phylum, sep = ' ')), size=5)+ labs(x="Phylum")+
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