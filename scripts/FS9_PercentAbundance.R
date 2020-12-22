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

#Reorder days -3 to 7 in 'fobar.gather' plot
fobar.gather.order$Day <- factor(fobar.gather.order$Day, levels=c("DNEG3", "D0", "D4", "D7"))
levels(sample_data(fobar.gather.order)$Day) #"DNEG3" "D0"    "D4"    "D7"  

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

###################################################### Q2 Day 7 Order NONINFnm vs INFnm #####################################################

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

#Remove DNEG3, D0, D4, NONINFnm
fobar.gather.order.2 <- fobar.gather.order %>% 
    select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Order", "value") %>% 
    filter(Day == "D7") %>% 
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

#Day 7 Order Figure: only Coriobacteriales, Mollicutes_RF39, and Verrucomicrobiales have % total community figures that match with DESeq2 data.
#Therefore, I will make a DESeq2-only figure that highlights all orders for Q2 D4, D7 that pass the log2-fold change > 0.25.
#Make figure that shows enriched in "INFfeed", "INFinject" or depleted in "INFfeed", "INFinject"
OrderFig_D7_2 <- fobar.gather.order.2 %>% filter(Day == "D7") %>% 
    select("group", "Day", "Pig", "Treatment", "Sample.type", "All", "Order", "value2") %>% 
    #filter(Order %in% c("Mollicutes_RF39", "Verrucomicrobiales", "Coriobacteriales")) %>% 
    ggplot(aes(x=Treatment, y=value2, group=All, fill=Treatment)) +
    geom_boxplot(position = 'identity') +
    geom_jitter(shape=21, width = .15)+
    facet_wrap("Order", scales = "free") +
    ylab('Percent of Total Community') +
    xlab ('') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x = element_blank()) +
          #axis.text.x=element_text(angle=45, hjust=1),
          #legend.text = element_text(face = "italic")) +
    guides(fill= guide_legend(ncol = 1)) +
    scale_fill_manual(values = c(INFnm='#CC0066', INFfeed='#999999', INFinject='#E69F00')) +
    theme_bw()
OrderFig_D7_2