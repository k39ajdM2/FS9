############################################################
#FS9 total phyla, orders
#Kathy Mou

#For generating total genera found in each treatment group per day and plot as bar graphs

#Clear workspace and load necessary packages
rm(list=ls())

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library("ggsci")
library(data.table)

############################################################

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


#################################################### Phylum #####################################################
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
    filter(Phylum != "Kiritimatiellaeota" & Phylum != "Verrucomicrobia") %>% 
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

###################################################### Order #####################################################
FS9.order <- tax_glom(FS9, 'Order')
order_tab <- as.data.frame(t(FS9.order@otu_table)) #Transpose 'FS9.order' by "otu_table"
head(order_tab)
FS9.order@tax_table[,4] #Fourth column is Order column
colnames(order_tab) <- FS9.order@tax_table[,4] #Replace column names in order_tab from Otuxxxx with Order names
order_tab2 <- order_tab/rowSums(order_tab) #Calculate the proportion of specific order per order column in 'order_tab'
head(order_tab2)
order_tab2$group <- rownames(order_tab2) #Create new column called "group" in 'order_tab2' containing rownames
head(order_tab2)
fobar.order <- merge(meta, order_tab2, by = 'group')
head(fobar.order)
fobar.gather.order <- fobar.order %>% gather(Order, value, -(group:All))
head(fobar.gather.order)

#Reorder days 0-14 in 'fobar.gather' plot
levels(sample_data(fobar.gather.order)$Day) #"D0"    "D4"    "D7"    "DNEG3"
fobar.gather.order$Day <- factor(fobar.gather.order$Day, levels=c("DNEG3", "D0", "D4", "D7"))
levels(sample_data(fobar.gather.order)$Day) #"DNEG3" "D0"    "D4"    "D7"  

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of order
fobar.gather.order %>% summarise_each(funs(n_distinct)) #48 total unique order
fobar.gather.order <- fobar.gather.order %>% group_by(All) %>% mutate(value2=(value/(length(All)/48))*100) %>% 
    arrange((desc(value2))) %>% 
    filter(value2 > 0) %>% 
    mutate_if(is.numeric, round, digits = 4)%>% 
    arrange(desc(value2)) %>% 
    ungroup()
    #48 refers to number of Order 

#Order Figures

#Day -3 Order
DNEG3Order <- fobar.gather.order %>% 
    subset(Day == "DNEG3") %>% 
    group_by(Order) %>% 
    select(Order, Treatment) %>% 
    count(Treatment)
#Decide which order to remove from plot that isn't present in all 4 treatment groups

OrderFig_DNEG3 <- fobar.gather.order %>% filter(Day == 'DNEG3') %>% 
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Order)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15)+
    facet_wrap(~Order, scales = 'free') + 
    ylab('Percent of Total Community') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Day -3") +
    scale_fill_igv(name = "Order") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 2))
OrderFig_DNEG3

#Day 0 Order
D0Order <- fobar.gather.order %>% 
    subset(Day == "D0") %>% 
    group_by(Order) %>% 
    select(Order, Treatment) %>% 
    count(Treatment)
#Decide which order to remove from plot that isn't present in all 4 treatment groups

OrderFig_D0 <- fobar.gather.order %>% filter(Day == 'D0') %>%
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Order)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15)+
    facet_wrap(~Order, scales = 'free') + 
    ylab('Percent of Total Community') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Day 0") +
    scale_fill_igv(name = "Order") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 2))
OrderFig_D0

#Day 4 Order
D4Order <- fobar.gather.order %>% 
    subset(Day == "D4") %>% 
    group_by(Order) %>% 
    select(Order, Treatment) %>% 
    count(Treatment)
#Decide which order to remove from plot that isn't present in all 4 treatment groups

OrderFig_D4 <- fobar.gather.order %>% filter(Day == 'D4') %>%
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Order)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15)+
    facet_wrap(~Order, scales = 'free') + 
    ylab('Percent of Total Community') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Day 4") +
    scale_fill_igv(name = "Order") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 2))
OrderFig_D4

#Day 7 Order
D7Order <- fobar.gather.order %>% 
    subset(Day == "D7") %>% 
    group_by(Order) %>% 
    select(Order, Treatment) %>% 
    count(Treatment)
#Decide which order to remove from plot that isn't present in all 4 treatment groups

OrderFig_D7 <- fobar.gather.order %>% filter(Day == 'D7' & value2 > 0) %>%
    filter(Order!= "Actinomycetales" & Order!="Bacillales" & Order!="Bacteroidetes_unclassified" & Order!="Betaproteobacteriales" &
               Order!= "Clostridia_unclassified" & Order!="Deltaproteobacteria_unclassified" & Order!="Elusimicrobiales" &
               Order!="Fibrobacterales" & Order!="Fusobacteriales" & Order!="Gammaproteobacteria_unclassified" & Order!="Pasteurellales" & Order!="Pyrinomonadales" & 
               Order!="Rhodospirillales" & Order!="Subgroup_6_or" & Order!="Synergistales" &
               Order!="Verrucomicrobiales" & Order!="WCHB1-41") %>% 
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Order)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15)+
    facet_wrap(~Order, scales = 'free') + 
    ylab('Percent of Total Community') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Day 7") +
    scale_fill_igv(name = "Order") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 1))
OrderFig_D7

write.csv(fobar.gather.order, file = "FS9_Order_OutDoubletons.csv")
