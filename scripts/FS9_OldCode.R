#Old code that I deleted from various R files


############################################################

#FS9_Percent Abundance.R

#FS9 percent abundance of total community - phyla level
#Kathy Mou

#For generating total genera found in each treatment group per day and plot as bar graphs

#Clear workspace and load necessary packages
rm(list=ls())

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("ggsci")
library(data.table)
library(cowplot)
library(wesanderson)

############################################################

#Must run this section first before running the other sections!

otu <- import_mothur(mothur_shared_file = './data/stability.outdoubletons.abund.opti_mcc.0.03.subsample.shared')
taxo <- import_mothur(mothur_constaxonomy_file = './data/stability.outdoubletons.abund.opti_mcc.0.03.cons.taxonomy')
meta <- read.table('./data/FS9_metadata.csv', header = TRUE, sep = ",")
meta$All <- with(meta, paste0(Day, sep="_", Treatment))
colnames(meta)[1] <- 'group' 
#Rename first column of "meta" as "group" temporarily. Will use "group" to set as rownames later and remove the "group" column
meta$group <- as.character(meta$group)
head(meta)
phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
head(phy_meta)
phy_meta <- phy_meta[,-1]
head(phy_meta)

#Create phyloseq-class objects with "otu" and "taxo"
FS9 <- phyloseq(otu, taxo)
FS9 <- merge_phyloseq(FS9, phy_meta)  #This combines the 'phy_meta' metadata with 'FS9' phyloseq object
colnames(tax_table(FS9)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
sample_sums(FS9) #Calculate the sum of all OTUs for each sample. All samples have 1300 sequences
FS9 <- prune_taxa(taxa_sums(FS9) > 2, FS9)  #Removes OTUs that occur less than 2 times globally

#################################################### Phylum All 4 groups #####################################################
FS9.phylum <- tax_glom(FS9, 'Phylum')
phyla_tab <- as.data.frame(t(FS9.phylum@otu_table)) #Transpose 'FS9.phylum' by "otu_table"
head(phyla_tab)
FS9.phylum@tax_table[,2] #Second column is Phylum column
colnames(phyla_tab) <- FS9.phylum@tax_table[,2] #Replace column names in phyla_tab from Otuxxxx with Phylum names
phyla_tab2 <- phyla_tab/rowSums(phyla_tab) #Calculate the proportion of specific phyla per phyla column in 'phyla_tab'
head(phyla_tab2)
phyla_tab3 <- phyla_tab2[,colSums(phyla_tab2)>0.1] #Keep the columns that have greater than 0.1 value
phyla_tab3$group <- rownames(phyla_tab3) #Rename rownames as "group"
head(phyla_tab3)
fobar <- merge(meta, phyla_tab3, by = 'group')
head(fobar)
fobar.gather <- fobar %>% gather(Phylum, value, -(group:All))  #Combines all phyla into one column and their respective % abundance values in "value" column
#It added columns "group" through "All" before "Phylum" and "value"
head(fobar.gather)

#Reorder days in 'fobar.gather' plot
levels(sample_data(fobar.gather)$Day) #"D0"    "D4"    "D7"    "DNEG3"
fobar.gather$Day <- factor(fobar.gather$Day, levels=c("DNEG3", "D0", "D4","D7"))
levels(sample_data(fobar.gather)$Day) #"DNEG3" "D0"    "D4"    "D7" 

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of phylum
fobar.gather %>% summarise_each(funs(n_distinct)) #13 total unique phyla
fobar.gather <- fobar.gather %>% 
  group_by(All) %>% 
  mutate(value2=(value/(length(All)/13))*100) %>% #13 refers to number of Phyla
  mutate(forplot=if_else(value2>0, "keep", "toss")) %>% 
  group_by(All) %>% 
  filter(forplot=='keep') %>% 
  mutate_if(is.numeric, round, digits = 4)%>% 
  arrange(desc(value2)) %>% 
  ungroup()

#Phylum Figures

#Day -3 Phylum
DNEG3Phylum <- fobar.gather %>% 
  subset(Day == "DNEG3") %>% 
  group_by(Phylum) %>% 
  select(Phylum, All) %>% 
  count(All)
#Decide which phyla to remove from plot that isn't present in all 4 treatment groups

PhylumFig_DNEG3 <- fobar.gather %>% filter(Day == 'DNEG3' & forplot == "keep") %>%
  ggplot(aes(x=Treatment, y=value2, group=All, fill=Phylum)) +
  geom_boxplot(position = 'identity') +
  geom_jitter(shape=21, width = .15)+
  facet_wrap(~Phylum, scales = 'free') + 
  ylab('Percent of Total Community') +
  xlab ('') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Day -3") +
  scale_fill_igv(name = "Phylum") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_text(size=5),
        legend.text = element_text(face = "italic"))
PhylumFig_DNEG3

#Day 0 Phylum
D0Phylum <- fobar.gather %>% 
  subset(Day == "D0") %>% 
  group_by(Phylum) %>% 
  select(Phylum, Treatment) %>% 
  count(Treatment)
#Decide which phyla to remove from plot that isn't present in all 4 treatment groups

PhylumFig_D0 <- fobar.gather %>% filter(Day == 'D0'  & forplot == "keep") %>%
  ggplot(aes(x=Treatment, y=value2, group=All, fill=Phylum)) +
  geom_boxplot(position = 'identity') +
  geom_jitter(shape=21, width = .15)+
  facet_wrap(~Phylum, scales = 'free') + 
  ylab('Percent of Total Community') +
  xlab ('') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Day 0") +
  scale_fill_igv(name = "Phylum") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"))
PhylumFig_D0

#Day 4 Phylum
D4Phylum <- fobar.gather %>% 
  subset(Day == "D4") %>% 
  group_by(Phylum) %>% 
  select(Phylum, Treatment) %>% 
  count(Treatment)
#Decide which phyla to remove from plot that isn't present in all 4 treatment groups

PhylumFig_D4 <- fobar.gather %>% filter(Day == 'D4' & forplot == "keep") %>%
  ggplot(aes(x=Treatment, y=value2, group=All, fill=Phylum)) +
  geom_boxplot(position = 'identity') +
  geom_jitter(shape=21, width = .15)+
  facet_wrap(~Phylum, scales = 'free') + 
  ylab('Percent of Total Community') +
  xlab ('') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Day 4") +
  scale_fill_igv(name = "Phylum") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"))
PhylumFig_D4

#Day 7 Phylum
D7Phylum <- fobar.gather %>% 
  subset(Day == "D7") %>% 
  group_by(Phylum) %>% 
  select(Phylum, Treatment) %>% 
  count(Treatment)
#Decide which phyla to remove from plot that isn't present in all 4 treatment groups

PhylumFig_D7 <- fobar.gather %>% filter(Day == 'D7' & forplot == "keep") %>%
  filter(Phylum != "Verrucomicrobia") %>% 
  ggplot(aes(x=Treatment, y=value2, group=All, fill=Phylum)) +
  geom_boxplot(position = 'identity') +
  geom_jitter(shape=21, width = .15)+
  facet_wrap(~Phylum, scales = 'free') + 
  ylab('Percent of Total Community') +
  xlab ('') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Day 7") +
  scale_fill_igv(name = "Phylum") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"))
PhylumFig_D7

write.csv(fobar.gather, file = "FS9_Phylum_OutDoubletons.csv")

#################################################### Q1 Days -3 and 7 Phylum NONINFnm vs INFnm #####################################################
FS9.phylum <- tax_glom(FS9, 'Phylum')
phyla_tab <- as.data.frame(t(FS9.phylum@otu_table)) #Transpose 'FS9.phylum' by "otu_table"
head(phyla_tab)
FS9.phylum@tax_table[,2] #Second column is Phylum column
colnames(phyla_tab) <- FS9.phylum@tax_table[,2] #Replace column names in phyla_tab from Otuxxxx with Phylum names
phyla_tab2 <- phyla_tab/rowSums(phyla_tab) #Calculate the proportion of specific phyla per phyla column in 'phyla_tab'
head(phyla_tab2)
phyla_tab3 <- phyla_tab2[,colSums(phyla_tab2)>0.1] #Keep the columns that have greater than 0.1 value
phyla_tab3$group <- rownames(phyla_tab3) #Rename rownames as "group"
head(phyla_tab3)
fobar <- merge(meta, phyla_tab3, by = 'group')
head(fobar)
fobar.gather <- fobar %>% gather(Phylum, value, -(group:All))  #Combines all phyla into one column and their respective % abundance values in "value" column
#It added columns "group" through "All" before "Phylum" and "value"
head(fobar.gather)

#Reorder days in 'fobar.gather' plot
levels(sample_data(fobar.gather)$Day) #"D0"    "D4"    "D7"    "DNEG3"
fobar.gather$Day <- factor(fobar.gather$Day, levels=c("DNEG3", "D0", "D4","D7"))
levels(sample_data(fobar.gather)$Day) #"DNEG3" "D0"    "D4"    "D7" 

#Remove INFinject, INFfeed
fobar.gather.phyla.q1 <- fobar.gather %>% 
  select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Phylum", "value") %>% 
  filter(Treatment != "INFfeed") %>% 
  filter(Treatment != "INFinject")

unique(fobar.gather.phyla.q1$Treatment) #"NONINFnm" "INFnm"

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of phylum
fobar.gather.phyla.q1 %>% summarise(n_distinct(fobar.gather.phyla.q1$Phylum)) #13 total unique phyla
fobar.gather.phyla.q1 <- fobar.gather.phyla.q1 %>% 
  group_by(All) %>% 
  mutate(value2=(value/(length(All)/13))*100) %>% #13 refers to number of Phyla
  arrange((desc(value2))) %>% 
  filter(value2 > 0) %>% 
  mutate_if(is.numeric, round, digits = 4)%>% 
  arrange(desc(value2)) %>% 
  ungroup()

#Phylum Figures
#Day -3 Phylum
DNEG3Phylum <- fobar.gather.phyla.q1 %>% 
  subset(Day == "DNEG3") %>% 
  group_by(Phylum) %>% 
  select(Phylum, All) %>% 
  count(All)
#Decide which phyla to remove from plot that isn't present in all 4 treatment groups

PhylumFig_DNEG3 <- fobar.gather.phyla.q1 %>% filter(Day == 'DNEG3') %>%
  select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Phylum", "value2") %>% 
  #    filter(Phylum %in% c("Actinobacteria", "Cyanobacteria", "Epsilonbacteraeota", "Firmicutes", "Spirochaetes")) %>% 
  ggplot(aes(x=Treatment, y=value2, group=All, fill=Treatment)) +
  geom_boxplot(position = 'identity') +
  geom_jitter(shape=21, width = .15)+
  facet_wrap(~Phylum, scales = 'free') + 
  ylab('Percent of Total Community (Phylum)') +
  xlab ('') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_text(size=5),
        legend.text = element_text(face = "italic")) +
  theme_bw()
PhylumFig_DNEG3

ggsave("Q1_NONINFnm_INFnm_Phylum_PercentAbundance_WithDeSeq2Data.tiff", plot=PhylumFig_DNEG3, width = 10, height = 6, dpi = 500, units =c("in"))

#Day 7 Phylum
D7Phylum <- fobar.gather.phyla.q1 %>% 
  subset(Day == "D7") %>% 
  group_by(Phylum) %>% 
  select(Phylum, Treatment) %>% 
  count(Treatment)
#Decide which phyla to remove from plot that isn't present in all 4 treatment groups

#Couldn't get NONINFnm to show up with INFnm for Verrucomicrobia (even when I set filter(value2 > 0)), so I will ignore PhylumFig_D7 plot
#PhylumFig_D7 <- fobar.gather.phyla.q1 %>% filter(Day == 'D7') %>%
#    select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Phylum", "value2") %>% 
#    filter(Phylum %in% c("Verrucomicrobia")) %>% 
#    ggplot(aes(x=Treatment, y=value2, group=All, fill=Treatment)) +
#    geom_boxplot(position = 'identity') +
#    geom_jitter(shape=21, width = .15)+
#    facet_wrap(~Phylum, scales = 'free') + 
#    ylab('Percent of Total Community') +
#    xlab ('') +
#    theme(plot.title = element_text(hjust = 0.5)) +
#    scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
#    theme(axis.text.x=element_text(angle=45, hjust=1),
#          axis.title.x = element_blank(),
#          legend.text = element_text(face = "italic")) +
#    theme_bw()
#PhylumFig_D7

#Re-running PhylumFig_D7 to see what it looks like if filter set to value2 > 0 vs value2 > 0.01
PhylumFig_D7 <- fobar.gather.phyla.q1 %>% filter(Day == 'D7') %>%
  select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Phylum", "value2") %>% 
  #filter(Phylum %in% c("Verrucomicrobia")) %>% 
  ggplot(aes(x=Treatment, y=value2, group=All, fill=Treatment)) +
  geom_boxplot(position = 'identity') +
  geom_jitter(shape=21, width = .15)+
  facet_wrap(~Phylum, scales = 'free') + 
  ylab('Percent of Total Community') +
  xlab ('') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9')) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic")) +
  theme_bw()
PhylumFig_D7


#Cowplot of Q1_NONINFnm_INFnm DNEG3 Order and Phylum
Fig2ab <- plot_grid(PhylumFig_DNEG3, OrderFig_DNEG3_Q1, labels = c('A', 'B'), label_size = 12, rel_widths = c(6, 5))
Fig2ab

ggsave("Q1_NONINFnm_INFnm_Phylum_Order_PercentAbundance_WithDeSeq2Data.tiff", plot=Fig2ab, width = 12, height = 8, dpi = 500, units =c("in"))




#####################################################################################################
#FS9 DESeq2 - Noninfected (NONINFnm) vs Infected (INFnm) only on all days, or days 4 and 7
#Kathy Mou

#Purpose: This code uses DESeq2 package to identify fecal microbial genera that were differentially 
#abundant between the two groups NONINFnm and INFnm on all days and only days 4 and 7 at the phylum and order level

#Load library packages
library(DESeq2)
library(phyloseq)
library(ggplot2)
library("wesanderson")
library(plotly)
library(gapminder)
library("ggsci")
library("apeglm")

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
FS9.phylum <- tax_glom(FS9, taxrank = "Phylum")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 


##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (PHYLUM LEVEL) BETWEEN TREATMENTS ###################################

##################################### Day -3, 0, 4, 7 ##########################################################################

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

######### Day -3 INFnm vs NONINFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "DNEG3_INFnm")
#20
sum(meta$Set == "DNEG3_NONINFnm")
#20

#re-level your factor and re-run DESeq2
sample_data(FS9.DNEG3.p)$Set <- factor(sample_data(FS9.DNEG3.p)$Set,
                                       levels =c('DNEG3_NONINFnm',"DNEG3_INFnm"))
FS9.DNEG3.p.De <- phyloseq_to_deseq2(FS9.DNEG3.p, ~ Set)
FS9.DNEG3.p.De <- DESeq(FS9.DNEG3.p.De, test = "Wald", fitType = "parametric")
FS9.DNEG3.p.De$Set
resultsNames(FS9.DNEG3.p.De)
#[1] "Intercept"                   "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm"  
res.DNEG3.p = lfcShrink(FS9.DNEG3.p.De, coef = "Set_DNEG3_INFnm_vs_DNEG3_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (DNEG3_INF)
# negative log2foldchanges are associated with the second group from this line (DNEG3_NONINF)
sigtab.DNEG3.p = res.DNEG3.p[which(res.DNEG3.p$padj < .05), ]
sigtab.DNEG3.p = cbind(as(sigtab.DNEG3.p, "data.frame"), as(tax_table(FS9.DNEG3.p)[rownames(sigtab.DNEG3.p), ], "matrix"))
sigtab.DNEG3.p$newp <- format(round(sigtab.DNEG3.p$padj, digits = 3), scientific = TRUE)
sigtab.DNEG3.p$Treatment <- ifelse(sigtab.DNEG3.p$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.DNEG3.p)

deseq.DNEG3.p <- 
  ggplot(sigtab.DNEG3.p, aes(x=reorder(rownames(sigtab.DNEG3.p), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.DNEG3.p), y=0, label = paste(Phylum)), size=5, fontface = 'italic')+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phylum\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day -3')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.DNEG3.p

#Add OTU and comparisons columns
sigtab.DNEG3.p
sigtab.DNEG3.p$OTU <- rownames(sigtab.DNEG3.p)
sigtab.DNEG3.p
sigtab.DNEG3.p$comp <- 'DNEG3_INFnm_vs_NONINFnm'

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

######### Day 0 INFinject vs INFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D0_INFnm")
#8
sum(meta$Set == "D0_NONINFnm")
#9

sample_data(FS9.D0.p)$Set <- factor(sample_data(FS9.D0.p)$Set,
                                    levels =c('D0_NONINFnm',"D0_INFnm"))
FS9.D0.p.De <- phyloseq_to_deseq2(FS9.D0.p, ~ Set)
FS9.D0.p.De <- DESeq(FS9.D0.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D0.p.De)
#[1] "Intercept"      "Set_D0_INFnm_vs_D0_NONINFnm"
res.D0.p = lfcShrink(FS9.D0.p.De, coef = "Set_D0_INFnm_vs_D0_NONINFnm", type = 'apeglm')
sigtab.D0.p = res.D0.p[which(res.D0.p$padj < .05), ]
sigtab.D0.p = cbind(as(sigtab.D0.p, "data.frame"), as(tax_table(FS9.D0.p)[rownames(sigtab.D0.p), ], "matrix"))
sigtab.D0.p$newp <- format(round(sigtab.D0.p$padj, digits = 3), scientific = TRUE)
sigtab.D0.p$Treatment <- ifelse(sigtab.D0.p$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D0.p) 

deseq.D0.p <- 
  ggplot(sigtab.D0.p, aes(x=reorder(rownames(sigtab.D0.p), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.p), y=0, label = paste(Phylum)), size=5, fontface = 'italic')+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phylum\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day 0')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.D0.p

#Add OTU and comparisons columns
sigtab.D0.p
sigtab.D0.p$OTU <- rownames(sigtab.D0.p)
sigtab.D0.p
sigtab.D0.p$comp <- 'D0_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab.phylum <- rbind(final.sigtab.phylum, sigtab.D0.p)




##################################################### Day 4 ######################################################################

FS9.D4.p <- subset_samples(FS9.phylum, Day == 'D4')
sample_sums(FS9.D4.p)
colnames(otu_table(FS9.D4.p)) #check on all the sample names
FS9.D4.p <- prune_taxa(taxa_sums(FS9.D4.p) > 1, FS9.D4.p)
#if taxa_sums is >1, then it will print that out in FS9.D4 object and not include anything with <1.
rowSums(FS9.D4.p@otu_table)

#Look at what Set is
unique(sample_data(FS9.D4.p)$Set)

######### Day 4 INFinject vs INFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D4_INFnm")
#6
sum(meta$Set == "D4_NONINFnm")
#5

sample_data(FS9.D4.p)$Set <- factor(sample_data(FS9.D4.p)$Set,
                                    levels =c('D4_NONINFnm',"D4_INFnm"))
FS9.D4.p.De <- phyloseq_to_deseq2(FS9.D4.p, ~ Set)
FS9.D4.p.De <- DESeq(FS9.D4.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D4.p.De)
#[1] "Intercept"      "Set_D4_INFnm_vs_D4_NONINFnm"
res.D4.p = lfcShrink(FS9.D4.p.De, coef = "Set_D4_INFnm_vs_D4_NONINFnm", type = 'apeglm')
sigtab.D4.p = res.D4.p[which(res.D4.p$padj < .05), ]
sigtab.D4.p = cbind(as(sigtab.D4.p, "data.frame"), as(tax_table(FS9.D4.p)[rownames(sigtab.D4.p), ], "matrix"))
sigtab.D4.p$newp <- format(round(sigtab.D4.p$padj, digits = 3), scientific = TRUE)
sigtab.D4.p$Treatment <- ifelse(sigtab.D4.p$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D4.p) 
#DataFrame with 0 rows and 7 columns. Will skip to next data set.




##################################################### Day 7 ######################################################################

FS9.D7.p <- subset_samples(FS9.phylum, Day == 'D7')
sample_sums(FS9.D7.p)
colnames(otu_table(FS9.D7.p)) #check on all the sample names
FS9.D7.p <- prune_taxa(taxa_sums(FS9.D7.p) > 1, FS9.D7.p)
#if taxa_sums is >1, then it will print that out in FS9.D7 object and not include anything with <1.
rowSums(FS9.D7.p@otu_table)

#Look at what Set is
unique(sample_data(FS9.D7.p)$Set)

######### Day 7 INFinject vs INFnm ###################

meta$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta$Set == "D7_INFnm")
#14
sum(meta$Set == "D7_NONINFnm")
#12

sample_data(FS9.D7.p)$Set <- factor(sample_data(FS9.D7.p)$Set,
                                    levels =c('D7_NONINFnm',"D7_INFnm"))
FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.p.De)
#[1] "Intercept"      "Set_D7_INFnm_vs_D7_NONINFnm"
res.D7.p = lfcShrink(FS9.D7.p.De, coef = "Set_D7_INFnm_vs_D7_NONINFnm", type = 'apeglm')
sigtab.D7.p = res.D7.p[which(res.D7.p$padj < .05), ]
sigtab.D7.p = cbind(as(sigtab.D7.p, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p), ], "matrix"))
sigtab.D7.p$newp <- format(round(sigtab.D7.p$padj, digits = 3), scientific = TRUE)
sigtab.D7.p$Treatment <- ifelse(sigtab.D7.p$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D7.p) 

deseq.D7.p <- 
  ggplot(sigtab.D7.p, aes(x=reorder(rownames(sigtab.D7.p), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p), y=0, label = paste(Phylum)), size=5, fontface = 'italic')+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phylum\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.D7.p

#Add OTU and comparisons columns
sigtab.D7.p
sigtab.D7.p$OTU <- rownames(sigtab.D7.p)
sigtab.D7.p
sigtab.D7.p$comp <- 'D7_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab.phylum <- rbind(final.sigtab.phylum, sigtab.D7.p)

#write csv
write.csv(final.sigtab.phylum, file= "FS9_FinalDiffAbund_Phylum_OutDoubletons_Q1_alldays.csv")




##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (PHYLUM LEVEL) BETWEEN TREATMENTS ###################################

##################################### Days 4 and 7 ##########################################################################

######################################################### Day 4 #########################################################

sample_data(FS9.phylum)

FS9.D4.p <- subset_samples(FS9.phylum, Day == 'D4')
sample_sums(FS9.D4.p)
colnames(otu_table(FS9.D4.p)) #check on all the sample names
FS9.D4.p <- prune_taxa(taxa_sums(FS9.D4.p) > 1, FS9.D4.p)
#if taxa_sums is >1, then it will print that out in FS9.D4.p object and not include anything with <1.
rowSums(FS9.D4.p@otu_table)

#Look at what Set is
sample_data(FS9.D4.p)

######### Day 4 INFnm vs NONINFnm ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D4_INFnm")
#6
sum(meta2$Set == "D4_NONINFnm")
#5

#re-level your factor and re-run DESeq2
sample_data(FS9.D4.p)$Set <- factor(sample_data(FS9.D4.p)$Set,
                                    levels =c('D4_NONINFnm',"D4_INFnm"))
FS9.D4.p.De <- phyloseq_to_deseq2(FS9.D4.p, ~ Set)
FS9.D4.p.De <- DESeq(FS9.D4.p.De, test = "Wald", fitType = "parametric")
FS9.D4.p.De$Set
resultsNames(FS9.D4.p.De)
#[1] "Intercept"                   "Set_D4_INFnm_vs_D4_NONINFnm"  
res.D4.p = lfcShrink(FS9.D4.p.De, coef = "Set_D4_INFnm_vs_D4_NONINFnm", type = 'apeglm')
# positive log2foldchanges are associated with the first group from this line (D4_INF)
# negative log2foldchanges are associated with the second group from this line (D4_NONINF)
sigtab.D4.p = res.D4.p[which(res.D4.p$padj < .05), ]
sigtab.D4.p = cbind(as(sigtab.D4.p, "data.frame"), as(tax_table(FS9.D4.p)[rownames(sigtab.D4.p), ], "matrix"))
sigtab.D4.p$newp <- format(round(sigtab.D4.p$padj, digits = 3), scientific = TRUE)
sigtab.D4.p$Treatment <- ifelse(sigtab.D4.p$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D4.p)
#DataFrame with 0 rows and 7 columns. Skip to next dataset.



##################################################### Day 7 ######################################################################

FS9.D7.p <- subset_samples(FS9.phylum, Day == 'D7')
sample_sums(FS9.D7.p)
colnames(otu_table(FS9.D7.p)) #check on all the sample names
FS9.D7.p <- prune_taxa(taxa_sums(FS9.D7.p) > 1, FS9.D7.p)
#if taxa_sums is >1, then it will print that out in FS9.D7 object and not include anything with <1.
rowSums(FS9.D7.p@otu_table)

#Look at what Set is
unique(sample_data(FS9.D7.p)$Set)

######### Day 7 INFnm vs NONINFnm ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D7_INFnm")
#14
sum(meta2$Set == "D7_NONINFnm")
#12

sample_data(FS9.D7.p)$Set <- factor(sample_data(FS9.D7.p)$Set,
                                    levels =c('D7_NONINFnm',"D7_INFnm"))
FS9.D7.p.De <- phyloseq_to_deseq2(FS9.D7.p, ~ Set)
FS9.D7.p.De <- DESeq(FS9.D7.p.De, test = "Wald", fitType = "parametric")
resultsNames(FS9.D7.p.De)
#[1] "Intercept"     "Set_D7_INFnm_vs_D7_NONINFnm"
res.D7.p = lfcShrink(FS9.D7.p.De, coef = "Set_D7_INFnm_vs_D7_NONINFnm", type = 'apeglm')
sigtab.D7.p = res.D7.p[which(res.D7.p$padj < .05), ]
sigtab.D7.p = cbind(as(sigtab.D7.p, "data.frame"), as(tax_table(FS9.D7.p)[rownames(sigtab.D7.p), ], "matrix"))
sigtab.D7.p$newp <- format(round(sigtab.D7.p$padj, digits = 3), scientific = TRUE)
sigtab.D7.p$Treatment <- ifelse(sigtab.D7.p$log2FoldChange >=0, "INFnm", "NONINFnm")
head(sigtab.D7.p) 

deseq.D7.p <- 
  ggplot(sigtab.D7.p, aes(x=reorder(rownames(sigtab.D7.p), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.p), y=0, label = paste(Phylum)), size=5, fontface = 'italic')+ labs(x="Phylum")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))+ ggtitle('Differentially Abundant Phylum\n in INFnm Group Relative to NONINFnm\n in Fecal Microbiota on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 14, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
deseq.D7.p

#Add OTU and comparisons columns
sigtab.D7.p
sigtab.D7.p$OTU <- rownames(sigtab.D7.p)
sigtab.D7.p
sigtab.D7.p$comp <- 'D7_INFnm_vs_NONINFnm'

#Create final significant comparisons table
final.sigtab.phylum <- sigtab.D7.p

#write csv
write.csv(final.sigtab.phylum, file= "FS9_FinalDiffAbund_Phylum_OutDoubletons_Q1_D4D7.csv")





#######################################################################################################
#FS9 DESeq2 - Infected: nm vs feed vs inject, days 4 and 7 - order level
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
FS9.phylum <- tax_glom(FS9, taxrank = "Phylum")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 

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

#write csv
write.csv(final.sigtab.phylum, file= "FS9_FinalDiffAbund_Phylum_OutDoubletons_Q2_D4D7.csv")


##################################################################################################################################

##################################### SIGNIFICANT CHANGES IN ABUNDANCE OF ORGANISMS (ORDER LEVEL) BETWEEN TREATMENTS ###################################

##################################################### Day 7 ######################################################################

######### 3. Day 7 INFfeed vs INFinject - took this out from FS9_DeSeq2_Q2.R because 
######### it was not relevant to the paper (paper will only look at INFinject or INFnm
########## relative to INFnm

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