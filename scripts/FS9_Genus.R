############################################################
#FS9 total phyla, orders
#Kathy Mou

#For generating total genera found in each treatment group per day and plot as bar graphs

#Set working directory
setwd("~/Desktop/FS9/FS9_RWorkspace")

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

#Load saved image
load("FS9.Genus.RData")



#Save image
save.image(file="FS9.Genus.RData")

############################################################

otu <- import_mothur(mothur_shared_file = 'stability.outsingletons.abund.opti_mcc.0.03.subsample.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'stability.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
shared <- read.table('stability.outsingletons.abund.opti_mcc.0.03.subsample.shared', header = TRUE)
meta <- read.table('FS9_metadata.csv', header = TRUE, sep = ",")
head(meta)

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
phyla_tab2$group <- rownames(phyla_tab2) #Create new column called "group" in 'phyla_tab2' containing rownames
head(phyla_tab2)

#Check to see if there is an extra "group" column. If so, run the next set of commands (up to "head(fobar)") and 
#remove appropriate column
which(colnames(phyla_tab2)=="group") #Results say column 20 is "group" column
phyla_tab3 <- phyla_tab2[,-20] #Drop the 20th column
phyla_tab4 <- phyla_tab3[,colSums(phyla_tab3)>0.1] #Keep the columns that have greater than 0.1 value
phyla_tab4$group <- rownames(phyla_tab4) #Rename rownames as "group"
fobar <- merge(meta, phyla_tab4, by = 'group')
head(fobar)

fobar.gather <- fobar %>% gather(Phylum, value, -(group:All))  #This converts 'fobar' to long-form dataframe. This is handy for using ggplot faceting functions, check out tidyverse tutorials
#This also created new columns "Phylum", "value"; it added columns "group" through "All" before "Phylum" and "value"
head(fobar.gather)

#Reorder days 0-14 in 'fobar.gather' plot
levels(sample_data(fobar.gather)$Day)
fobar.gather$Day <- factor(fobar.gather$Day, levels=c("-3", "0", "7"))
head(fobar.gather$Day)

#Create "All" column with "Day", "Treatment" and "Tissue" in 'fobar.gather'
fobar.gather$All <- paste(fobar.gather$Day, fobar.gather$Treatment, sep = '_')
head(fobar.gather$All)

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of phylum
fobar.gather %>% summarise_each(funs(n_distinct)) #19 total unique phyla
fobar.gather <- fobar.gather %>% group_by(All) %>% mutate(value2=(value/(length(All)/19))*100) #19 refers to number of Phyla

#Phylum Figures

#Day -3 Phylum
PhylumFig_DNEG3 <- fobar.gather %>% filter(Day == '-3') %>%
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
          axis.title.x = element_blank())
PhylumFig_DNEG3
ggsave("FS9_Phylum_DNEG3.tiff", plot=PhylumFig_DNEG3, width = 15, height = 10, dpi = 500, units =c("in"))

#Day 0 Phylum
PhylumFig_D0 <- fobar.gather %>% filter(Day == '0') %>%
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
          axis.title.x = element_blank())
PhylumFig_D0
ggsave("FS9_Phylum_D0.tiff", plot=PhylumFig_D0, width = 15, height = 10, dpi = 500, units =c("in"))

#Day 7 Phylum
PhylumFig_D7 <- fobar.gather %>% filter(Day == '7') %>%
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
          axis.title.x = element_blank())
PhylumFig_D7
ggsave("FS9_Phylum_D7.tiff", plot=PhylumFig_D7, width = 15, height = 10, dpi = 500, units =c("in"))



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

#Check to see if there is an extra "group" column. If so, run the next set of commands (up to "head(fobar)") and 
#remove appropriate column
which(colnames(order_tab2)=="group") #Results say column 50 is "group" column
order_tab3 <- order_tab2[,-50] #Drop the 50th column
order_tab4 <- order_tab3[,colSums(order_tab3)>0.1] #Keep the columns that have greater than 0.1 value
order_tab4$group <- rownames(order_tab4) #Rename rownames as "group"
fobar2 <- merge(meta, order_tab4, by = 'group')
head(fobar2)

fobar2.gather <- fobar2 %>% gather(Order, value, -(group:All))  #This converts 'fobar' to long-form dataframe. This is handy for using ggplot faceting functions, check out tidyverse tutorials
#This also created new columns "Order", "value"; it added columns "group" through "All" before "Order" and "value"
head(fobar2.gather)

#Reorder days 0-14 in 'fobar.gather' plot
levels(sample_data(fobar2.gather)$Day)
fobar2.gather$Day <- factor(fobar2.gather$Day, levels=c("-3", "0", "7"))
head(fobar2.gather$Day)

#Create "All" column with "Day", "Treatment" and "Tissue" in 'fobar.gather'
fobar2.gather$All <- paste(fobar2.gather$Day, fobar2.gather$Treatment, sep = '_')
head(fobar2.gather$All)

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of order
fobar2.gather %>% summarise_each(funs(n_distinct)) #22 total unique order
fobar2.gather <- fobar2.gather %>% group_by(All) %>% mutate(value2=(value/(length(All)/22))*100) #22 refers to number of Order

#Order Figures

#Day -3 Order
OrderFig_DNEG3 <- fobar2.gather %>% filter(Day == '-3') %>%
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
          axis.title.x = element_blank())
OrderFig_DNEG3
ggsave("FS9_Order_DNEG3.tiff", plot=OrderFig_DNEG3, width = 17, height = 10, dpi = 500, units =c("in"))

#Day 0 Order
OrderFig_D0 <- fobar2.gather %>% filter(Day == '0') %>%
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
          axis.title.x = element_blank())
OrderFig_D0
ggsave("FS9_Order_D0.tiff", plot=OrderFig_D0, width = 17, height = 10, dpi = 500, units =c("in"))

#Day 7 Order
OrderFig_D7 <- fobar2.gather %>% filter(Day == '7') %>%
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
          axis.title.x = element_blank())
OrderFig_D7
ggsave("FS9_Order_D7.tiff", plot=OrderFig_D7, width = 17, height = 10, dpi = 500, units =c("in"))

write.csv(fobar2.gather, file = "FS9_Order.csv")
