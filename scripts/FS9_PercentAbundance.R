############################################################
#FS9 percent abundance of total community - order level
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




###################################################### Order All 4 groups #####################################################
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

#Reorder days 4 to 14 in 'fobar.gather' plot
fobar.gather.order$Day <- factor(fobar.gather.order$Day, levels=c("D4", "D7", "D11", "D14"))
levels(sample_data(fobar.gather.order)$Day) #"D4"  "D7"  "D11" "D14"

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of order
fobar.gather.order %>% summarise(n_distinct(fobar.gather.order$Order)) #48 total unique order
fobar.gather.order <- fobar.gather.order %>% group_by(All) %>% mutate(value2=(value/(length(All)/48))*100) %>% 
    arrange((desc(value2))) %>% 
    filter(value2 > 0) %>% 
    mutate_if(is.numeric, round, digits = 4)%>% 
    arrange(desc(value2)) %>% 
    ungroup()
    #48 refers to number of Order 

#Order Figures

#Day 4 Order
D4Order <- fobar.gather.order %>% 
    subset(Day == "D4") %>% 
    group_by(Order) %>% 
    select(Order, Treatment) %>% 
    count(Treatment)

#Decide which order to remove from plot that isn't present in all 4 treatment groups

(OrderFig_D4 <- fobar.gather.order %>% filter(Day == 'D4') %>% 
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
    guides(fill= guide_legend(ncol = 2)))

#Day 7 Order
D7Order <- fobar.gather.order %>% 
    subset(Day == "D7") %>% 
    group_by(Order) %>% 
    select(Order, Treatment) %>% 
    count(Treatment)
#Decide which order to remove from plot that isn't present in all 4 treatment groups
write.csv(D7Order, file = "FS9_D7Order_OutDoubletons.csv")

(OrderFig_D7 <- fobar.gather.order %>% filter(Day == 'D7') %>%
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
    guides(fill= guide_legend(ncol = 2)))

#Day 11 Order
D11Order <- fobar.gather.order %>% 
    subset(Day == "D11") %>% 
    group_by(Order) %>% 
    select(Order, Treatment) %>% 
    count(Treatment)
#Decide which order to remove from plot that isn't present in all 4 treatment groups
write.csv(D11Order, file = "FS9_D11Order_OutDoubletons.csv")

(OrderFig_D11 <- fobar.gather.order %>% filter(Day == 'D11') %>%
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Order)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15)+
    facet_wrap(~Order, scales = 'free') + 
    ylab('Percent of Total Community') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Day 11") +
    scale_fill_igv(name = "Order") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 2)))

#Day 14 Order
D14Order <- fobar.gather.order %>% 
    subset(Day == "D14") %>% 
    group_by(Order) %>% 
    select(Order, Treatment) %>% 
    count(Treatment)
#Decide which order to remove from plot that isn't present in all 4 treatment groups
write.csv(D14Order, file = "FS9_D14Order_OutDoubletons.csv")

(OrderFig_D14 <- fobar.gather.order %>% filter(Day == 'D14' & value2 > 0) %>%
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Order)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15)+
    facet_wrap(~Order, scales = 'free') + 
    ylab('Percent of Total Community') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Day 14") +
    scale_fill_igv(name = "Order") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 1)))

write.csv(fobar.gather.order, file = "FS9_Order_OutDoubletons.csv")


###################################################### Q1 Day -3 Order NONINFnm vs INFnm #####################################################

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

#Reorder days -3 to 7 in 'fobar.gather' plot
fobar.gather.order$Day <- factor(fobar.gather.order$Day, levels=c("DNEG3", "D0", "D4", "D7"))
levels(sample_data(fobar.gather.order)$Day) #"DNEG3" "D0"    "D4"    "D7"  

#Remove INFfeed, INFinject
fobar.gather.order.q1 <- fobar.gather.order %>% 
    select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Order", "value") %>% 
    filter(Treatment != "INFfeed") %>% 
    filter(Treatment != "INFinject")

unique(fobar.gather.order.q1$Treatment) #"NONINFnm" "INFnm"

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of phylum
fobar.gather.order.q1 %>% summarise(n_distinct(fobar.gather.order.q1$Order)) #48 total unique order
fobar.gather.order.q1 <- fobar.gather.order.q1 %>% 
    group_by(All) %>% 
    mutate(value2=(value/(length(All)/48))*100) %>% #48 refers to number of order
    arrange((desc(value2))) %>% 
    filter(value2 > 0) %>% 
    mutate_if(is.numeric, round, digits = 4)%>% 
    arrange(desc(value2)) %>% 
    ungroup()

unique(fobar.gather.order.q1$Order) #43 unique orders

#Day -3 Order Figure
OrderFig_DNEG3_Q1 <- fobar.gather.order.q1 %>% filter(Day == "DNEG3") %>% 
    select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Order", "value2") %>% 
    filter(Order %in% c("Betaproteobacteriales", "Campylobacterales", "Gastranaerophilales", "Lactobacillales", "Pasteurellales")) %>% 
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Treatment)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15) +
    facet_wrap("Order", scales = "free") +
    ylab('Percent of Total Community (Order)') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(strip.text = element_text(size=17),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          legend.title=element_text(size=14),
          legend.text = element_text(size=14)) +
    guides(fill= guide_legend(ncol = 1)) +
    theme_bw() +
    scale_fill_manual(values = c(INFnm='#CC0066', NONINFnm='#56B4E9'))
OrderFig_DNEG3_Q1

ggsave("Q1_NONINFnm_INFnm_Order_DNEG3_PercentAbundance_WithDeSeq2Data.tiff", plot=OrderFig_DNEG3_Q1, width = 9, height = 6, dpi = 500, units =c("in"))

###################################################### Q2 Day 14 Order NONINFnm vs INFnm #####################################################

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

#Reorder days 4 to 14 in 'fobar.gather' plot
fobar.gather.order$Day <- factor(fobar.gather.order$Day, levels=c("D4", "D7", "D11", "D14"))
levels(sample_data(fobar.gather.order)$Day) #"D4"  "D7"  "D11" "D14"

#Remove D4, D7, D11, NONINFnm
fobar.gather.order.2 <- fobar.gather.order %>% 
    select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Order", "value") %>% 
    filter(Day == "D14") %>% 
    filter(Treatment != "NONINFnm")

unique(fobar.gather.order.2$Treatment) #"INFfeed"   "INFinject" "INFnm"

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of order
fobar.gather.order.2 %>% summarise(n_distinct(fobar.gather.order.2$Order)) #48 total unique order
fobar.gather.order.2 <- fobar.gather.order.2 %>% group_by(All) %>% mutate(value2=(value/(length(All)/48))*100) %>% 
    arrange((desc(value2))) %>% 
    filter(value2 > 0) %>% 
    mutate_if(is.numeric, round, digits = 4)%>% 
    arrange(desc(value2)) %>% 
    ungroup()
#48 refers to number of Order 

unique(fobar.gather.order.2$Order) #33 unique orders

#Day 14 Order Figure: only Coriobacteriales, Mollicutes_RF39, and Verrucomicrobiales have % total community figures that match with DESeq2 data.
#Therefore, I will make a DESeq2-only figure that highlights all orders for Q2 D4, D7 that pass the log2-fold change > 0.25.
#Make figure that shows enriched in "INFfeed", "INFinject" or depleted in "INFfeed", "INFinject"
(OrderFig_D14_2 <- fobar.gather.order.2 %>% filter(Day == "D14") %>% 
    select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Order", "value2") %>% 
    filter(Order %in% c("Mollicutes_RF39", "Verrucomicrobiales", "Coriobacteriales", "Gammaproteobacteria_unclassified", "Rhodospirillales")) %>% 
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Treatment)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    facet_wrap("Order", scales = "free") +
    ylab('Percent of Total Community') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(fill= guide_legend(ncol = 1)) +
    scale_fill_manual(values = c(INFinject='#00BA38',
                                 INFnm='#F8766D',
                                 INFfeed='#619CFF')) +
    theme_bw())


###################################################### Phylum 3 Infected groups #####################################################
FS9.phylum <- tax_glom(FS9, 'Phylum')
phylum_tab <- as.data.frame(t(FS9.phylum@otu_table)) #Transpose 'FS9.phylum' by "otu_table"
head(phylum_tab)
FS9.phylum@tax_table[,2] #Second column is Phylum column
colnames(phylum_tab) <- FS9.phylum@tax_table[,2] #Replace column names in phylum_tab from Otuxxxx with Phylum names
phylum_tab2 <- phylum_tab/rowSums(phylum_tab) #Calculate the proportion of specific Phylum per Phylum column in 'phylum_tab'
head(phylum_tab2)
phylum_tab2$group <- rownames(phylum_tab2) #Create new column called "group" in 'phylum_tab2' containing rownames
head(phylum_tab2)
fobar.phylum <- merge(meta, phylum_tab2, by = 'group')
head(fobar.phylum)
fobar.gather.phylum <- fobar.phylum %>% gather(Phylum, value, -(group:All))
head(fobar.gather.phylum)

#Reorder days -3 to 7 in 'fobar.gather' plot
fobar.gather.phylum$Day <- factor(fobar.gather.phylum$Day, levels=c("DNEG3", "D0", "D4", "D7"))
levels(sample_data(fobar.gather.phylum)$Day) #"DNEG3" "D0"    "D4"    "D7"  

#Count the number of unique items in 'fobar.gather'. We're interested in the total unique number of order
fobar.gather.phylum %>% summarise(n_distinct(fobar.gather.phylum$Phylum)) #19 total unique order
fobar.gather.phylum <- fobar.gather.phylum %>% group_by(All) %>% mutate(value2=(value/(length(All)/19))*100) %>% 
    arrange((desc(value2))) %>% 
    filter(value2 > 0) %>% 
    mutate_if(is.numeric, round, digits = 4)%>% 
    arrange(desc(value2)) %>% 
    ungroup()
#19 refers to number of Phylum 

#Phylum Figures

#Day 4 Phylum
D4Phylum <- fobar.gather.phylum %>% 
    subset(Day == "D4") %>% 
    group_by(Phylum) %>% 
    select(Phylum, Treatment) %>% 
    count(Treatment)
#Decide which Phylum to remove from plot that isn't present in all 4 treatment groups

PhylumFig_D4 <- fobar.gather.phylum %>% filter(Day == 'D4') %>%
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
          legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 2))
PhylumFig_D4

#Day 7 Phylum
D7Phylum <- fobar.gather.phylum %>% 
    subset(Day == "D7") %>% 
    group_by(Phylum) %>% 
    select(Phylum, Treatment) %>% 
    count(Treatment)
#Decide which Phylum to remove from plot that isn't present in all 4 treatment groups

PhylumFig_D7 <- fobar.gather.phylum %>% filter(Day == 'D7' & value2 > 0) %>%
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
          legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 1))
PhylumFig_D7


PhylumFig_Bacteroidetes <- fobar.gather.phylum %>% select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Phylum", "value2") %>% 
    filter(Day %in% c('D4', 'D7') & Phylum == "Bacteroidetes" & Treatment != "NONINFnm") %>% 
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Phylum)) +
    geom_boxplot(position = 'identity', show.legend = FALSE) +
    geom_jitter(shape=21, width = .15, show.legend = FALSE) +
    facet_wrap(~Day, scales = 'fixed') + 
    ylab('Percent of Total Community') +
    scale_y_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7, 8)) +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Bacteroidetes") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          plot.title = element_text(face="italic"))
PhylumFig_Bacteroidetes

PhylumFig_Firmicutes <- fobar.gather.phylum %>% select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Phylum", "value2") %>% 
    filter(Day %in% c('D4', 'D7') & Phylum == "Firmicutes" & Treatment != "NONINFnm") %>% 
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Phylum)) +
    geom_boxplot(position = 'identity', fill = "lightblue", show.legend = FALSE) +
    geom_jitter(shape=21, width = .15, fill="lightblue", show.legend = FALSE) +
    facet_wrap(~Day, scales = 'fixed') + 
    ylab('Percent of Total Community') +
    scale_y_continuous(breaks=c(5, 7.5, 10, 12.5, 15, 17.5, 20)) +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Firmicutes") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          plot.title = element_text(face="italic"))
PhylumFig_Firmicutes

PhylumFig_BF <- fobar.gather.phylum %>% select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Phylum", "value2") %>% 
    filter(Day %in% c('D4', 'D7') & Treatment != "NONINFnm") %>% 
    filter(Phylum %in% c('Bacteroidetes', 'Firmicutes')) %>% 
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Phylum)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15) +
    facet_wrap(Phylum~Day, scales = 'fixed') + 
    ylab('Percent of Total Community') +
    scale_y_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)) +
    xlab ('') +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          strip.text = element_text(face="italic"),
          legend.text = element_text(face="italic"))
PhylumFig_BF

ggsave("Q2_INFnm_feed_inject_Bacteroidetes_PercentAbundance.tiff", plot=PhylumFig_Bacteroidetes, width = 6, height = 4, dpi = 500, units =c("in"))
ggsave("Q2_INFnm_feed_inject_Firmicutes_PercentAbundance.tiff", plot=PhylumFig_Firmicutes, width = 6, height = 4, dpi = 500, units =c("in"))
ggsave("Q2_INFnm_feed_inject_Bacteroidetes_Firmicutes_PercentAbundance.tiff", plot=PhylumFig_BF, width = 5, height = 6, dpi = 500, units =c("in"))
