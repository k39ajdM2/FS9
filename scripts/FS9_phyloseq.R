#########################################
#FS9 16S phyloseq 
#Kathy Mou

#NOTES: 
#This code creates phyloseq objects ("All" column) for beta diversity measures

#Clear workspace and load necessary packages
rm(list=ls())

#Load library packages
library(vegan)
library(tidyverse)
library(phyloseq)

sessionInfo()
#R version 3.6.3 (2020-02-29)

#########################################

#Read files for metadata and OTU table
meta <- read.csv("./data/FS9_metadata.csv", row.names = 1)
otu <- read.csv("./data/FS9.OTUtable.csv", row.names=1)
dim(otu) #1583 169
head(otu[,165:169])

#Remove taxonomy from 'otu'
tax <- otu[,(168:169)] #Removed column 169 "Taxonomy" and 168 to include the row names
head(tax)
colnames(tax)[1] <- "delete" #Renamed column 1 (formerly 168) as "delete" which will later be deleted
head(tax)

#Modify 'otu' with only OTU count data
otu <- otu[,-169] #Remove column 169 "Taxonomy" to have only OTU data
head(otu[,165:168]) 
dim(otu) #1583 rows 168 columns

#Transpose 'otu' to match format of 'meta'
otu.trans <- t(otu) #Now rownames are sample names, columns are OTUs
head(otu.trans[,1:5])
class(meta) #dataframe
class(otu) #dataframe

#Merge 'otu' and 'meta' data frames
otu.meta <- merge(meta, otu.trans, by.x=0, by.y=0) 
#x=0 means match via rownames from 'meta'; y=0 means match via rownames from 'otu.trans'
head(otu.meta[,1:10])
dim(otu.meta) #166 1588
rownames(otu.meta) <- otu.meta[,1] #Set column 1 as rownames for 'otu.meta'
otu.meta <- otu.meta[,-1] #Remove first column "rownames"
class(otu.meta) #dataframe
rownames(otu.meta)

#Added an "All" column (combines 'Treatment' and 'Day' values) in new 'otu.meta2'
otu.meta2<- cbind(otu.meta) #Make second copy of otu.meta to use to include "All" column
colnames(otu.meta2)
otu.meta2$All <- with(otu.meta2, paste0(Day, sep="_", Treatment)) #Create "All" column with Day and Treatment combined
head(otu.meta2)
dim(otu.meta2) #166 1588
head(otu.meta2[,1584:1588])
head(otu.meta2[,1:10])
otu.meta2<- otu.meta2[,c(1:4,1588,5:1587)] #Reorder columns to have "All" column after "Treatment" column
write.csv(otu.meta2, file="FS9.otu.meta_all.csv")

#Pull out metadata from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
meta2 <- otu.meta2[,c(1:5)] #Take columns 1-5 ("Day" to "All") from 'otu.meta2' to make 'meta2'
head(meta2)
dim(meta2) #166 5

#Create SAM metadata table phyloseq object
SAM = sample_data(meta2, errorIfNULL = TRUE)
head(SAM)
dim(SAM) #166 5

#Pull out OTU data from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
tail(otu.meta2[,1584:1588])
otu2 <- otu.meta2[,c(6:1588)] #Select OTU columns to create 'otu2' dataframe
head(otu2[,1:10])
dim(otu2) #166 1583
otu2.trans <- t(otu2) #Transpose otu.all to have OTUs as rownames, sample names as column names
head(otu2.trans[,1:10])
dim(otu2.trans) #1583 166

#Merge 'tax' back into 'otu2.trans' for correct format and taxons
head(tax)
otu2.tax <- merge(otu2.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu2.tax) #1583 169
head(otu2.tax[,1:10])
head(otu2.tax[,160:169])
row.names(otu2.tax) <- otu2.tax[,1] #Set first row of 'otu2.tax' as rownames
head(otu2.tax[,1:5])
otu2.tax <- otu2.tax[,-1] #Remove first row aka extraneous "Row.names" column from 'otu2.tax'
head(otu2.tax[,1:5])

#Split again
dim(otu2.tax) #1583 168
head(otu2.tax[,160:168])
otu2.notax <- otu2.tax[,1:166] #take rows 1-166 to make new dataframe 'otu2.notax' (167 is "delete" column, 168 is "Taxonomy" column)
dim(otu2.notax) #1583 166
head(otu2.notax[,1:5])
head(otu2.notax[,160:166])
class(otu2.notax) #dataframe
otu2.notax <- as.matrix(otu2.notax) #turn 'otu2.notax' into a matrix class
class(otu2.notax) #matrix
otu2.notax
otu2.notax.trans <- t(otu2.notax)
head(otu2.notax.trans[,1:10])

#Create OTU table phyloseq object
OTU = otu_table(otu2.notax.trans, taxa_are_rows = FALSE)
head(OTU[,1:10])
dim(OTU) #166 1583
class(OTU)
OTU2 <- prune_taxa(taxa_sums(OTU) > 0, OTU)
class(OTU2)
head(OTU2[,1:10])
taxa_sums(OTU2)
dim(OTU2) #166 1582

#Edit taxonomy
dim(otu2.tax) #1583 168
head(otu2.tax[,160:168])
tax2 <- separate(data = otu2.tax, 
                           col = Taxonomy, 
                           into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#"separate" function separates "Taxonomy" column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax2) #notice that Species column is blank
dim(tax2) #1583 174
head(tax2[,168:173])
tax2.kg <- tax2[,168:173] #Keep only taxonomy columns "Kingdom" to "Genus"
head(tax2.kg)
dim(tax2.kg) #1583 6
class(tax2.kg) #dataframe
tax2.kg <- as.matrix(tax2.kg)
class(tax2.kg) #matrix
tax2.kg

#Create TAX taxonomy table phyloseq object
TAX = tax_table(tax2.kg)
head(TAX)
dim(TAX) #1583 6

#Create phyloseq object containing taxonomy, metadata, and otu table
phyloseq.FS9 <- phyloseq(OTU2, SAM, TAX) #prune out any OTUs that have total to 0 in all samples
phyloseq.FS9
#otu_table()   OTU Table:         [ 1582 taxa and 166 samples ]
#sample_data() Sample Data:       [ 166 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 1582 taxa by 6 taxonomic ranks ]

save(phyloseq.FS9, file="phyloseq.FS9.RData") #Use for FS9_alpha_beta_diversity.R when loading phyloseq.FS9 object

###############################################################################

#Distance calculation
phyloseq.vegdist <- vegdist(phyloseq.FS9@otu_table, method="bray") 
phyloseq.adonis <- as(sample_data(phyloseq.FS9), "data.frame")
set.seed(1)
adonis.FS9 <- adonis(phyloseq.vegdist~Day*Treatment, data=phyloseq.adonis, permutations=9999)
adonis.FS9

#Call:
#adonis(formula = phyloseq.vegdist ~ Day * Treatment, data = phyloseq.adonis,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                 Df SumsOfSqs  MeanSqs   F.Model   R2        Pr(>F)    
#Day               3     5.582  1.86051   9.6752    0.14888   0.0001 ***
#Treatment         3     1.246  0.41542   2.1603    0.03324   0.0001 ***
#Day:Treatment     9     1.819  0.20209   1.0509    0.04851   0.3053    
#Residuals      150     28.845  0.19230             0.76937           
#Total          165     37.491                      1.00000          

#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Only Day, Treatment had significant effects on the gut microbial community variation observed between four groups