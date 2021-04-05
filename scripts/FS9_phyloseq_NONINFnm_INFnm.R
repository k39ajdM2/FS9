#########################################
#FS9 16S phyloseq - Noninfected (NONINFnm) vs infected (INFnm) days 7, 11, 14
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

#Read files for metadata and OTU table
meta <- read.csv("./data/FS9_metadata_NONINFnm_vs_INFnm.csv", row.names = 1)
otu <- read.csv("./data/FS9.OTUtable.doubleton.csv", row.names=1)
dim(otu) #1405 168
head(otu[,165:168])

#Remove taxonomy from 'otu'
tax <- otu[,(167:168)] #Removed column 168 "Taxonomy" and 167 to include the row names
head(tax)
colnames(tax)[1] <- "delete" #Renamed column 1 (formerly 167) as "delete" which will later be deleted
head(tax)

#Modify 'otu' with only OTU count data
otu <- otu[,-168] #Remove column 168 "Taxonomy" to have only OTU data
head(otu[,165:167]) 
dim(otu) #1405 rows 167 columns

#Transpose 'otu' to match format of 'meta'
otu.trans <- t(otu) #Now rownames are sample names, columns are OTUs
head(otu.trans[,1:5])
class(meta) #dataframe
class(otu) #dataframe

#Merge 'otu' and 'meta' data frames
otu.meta <- merge(meta, otu.trans, by.x=0, by.y=0) 
#x=0 means match via rownames from 'meta'; y=0 means match via rownames from 'otu.trans'
head(otu.meta[,1:10])
dim(otu.meta) #45 1410 (doubletons removed)
rownames(otu.meta) <- otu.meta[,1] #Set column 1 as rownames for 'otu.meta'
otu.meta <- otu.meta[,-1]
class(otu.meta) #dataframe
rownames(otu.meta)

#Added an "All" column (combines 'Treatment' and 'Day' values) in new 'otu.meta2'
otu.meta2<- cbind(otu.meta) #Make second copy of otu.meta to use to include "All" column
colnames(otu.meta2)
otu.meta2$All <- with(otu.meta2, paste0(Day, sep="_", Treatment)) #Create "All" column with Day and Treatment combined
head(otu.meta2)
dim(otu.meta2) #45 1410 (doubletons removed)
head(otu.meta2[,1405:1410])
head(otu.meta2[,1:10])
otu.meta2<- otu.meta2[,c(1:4,1410,5:1409)] #Reorder columns to have "All" column after "Treatment" column
#write.csv(otu.meta2, file="FS9.otu.meta_all.doubleton.csv")

#Pull out metadata from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
meta2 <- otu.meta2[,c(1:5)] #Take columns 1-5 ("Day" to "All") from 'otu.meta2' to make 'meta2'
head(meta2)
dim(meta2) #45 5 (doubletons removed)

#Create SAM metadata table phyloseq object
SAM = sample_data(meta2, errorIfNULL = TRUE)
head(SAM)
dim(SAM) #45 5

#Pull out OTU data from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
tail(otu.meta2[,1405:1410])
otu2 <- otu.meta2[,c(6:1410)] #Select OTU columns to create 'otu2' dataframe
head(otu2[,1:10])
dim(otu2) #45 1405 (doubletons removed)
otu2.trans <- t(otu2) #Transpose otu.all to have OTUs as rownames, sample names as column names
head(otu2.trans[,1:10])
dim(otu2.trans) #1405 45 (doubletons removed)

#Merge 'tax' back into 'otu2.trans' for correct format and taxons
head(tax)
otu2.tax <- merge(otu2.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu2.tax) #1405 48 (doubletons removed)
head(otu2.tax[,1:10])
head(otu2.tax[,40:48])
row.names(otu2.tax) <- otu2.tax[,1] #Set first row of 'otu2.tax' as rownames
head(otu2.tax[,1:5])
otu2.tax <- otu2.tax[,-1] #Remove first row aka extraneous "Row.names" column from 'otu2.tax'
head(otu2.tax[,1:5])

#Split again
dim(otu2.tax) #1405 47 (doubletons removed)
head(otu2.tax[,40:47])
otu2.notax <- otu2.tax[,1:45] #take rows 1-29 to make new dataframe 'otu2.notax' (46 is "delete" column, 47 is "Taxonomy" column)
dim(otu2.notax) #1405 45 (doubletons removed)
head(otu2.notax[,1:5])
head(otu2.notax[,40:45])
class(otu2.notax) #dataframe
otu2.notax <- as.matrix(otu2.notax) #turn 'otu2.notax' into a matrix class
class(otu2.notax) #matrix
otu2.notax
otu2.notax.trans <- t(otu2.notax)
head(otu2.notax.trans[,1:10])

#Create OTU table phyloseq object
OTU = otu_table(otu2.notax.trans, taxa_are_rows = FALSE)
head(OTU[,1:10])
dim(OTU) #45 1405 (doubletons removed)
class(OTU)
OTU2 <- prune_taxa(taxa_sums(OTU) > 0, OTU)
class(OTU2)
head(OTU2[,1:10])
taxa_sums(OTU2)
dim(OTU2) #45 901 (doubletons removed)

#Edit taxonomy
dim(otu2.tax) #1405 47 (doubletons removed)
head(otu2.tax[,40:47])
tax2 <- separate(data = otu2.tax, 
                 col = Taxonomy, 
                 into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#"separate" function separates "Taxonomy" column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax2) #notice that Species column is blank
dim(tax2) #1405 53 (doubletons removed)
head(tax2[,45:53])
tax2.kg <- tax2[,47:52] #Keep only taxonomy columns "Kingdom" to "Genus"
head(tax2.kg)
dim(tax2.kg) #1405 6 (doubletons removed)
class(tax2.kg) #dataframe
tax2.kg <- as.matrix(tax2.kg)
class(tax2.kg) #matrix
tax2.kg

#Create TAX taxonomy table phyloseq object
TAX = tax_table(tax2.kg)
head(TAX)
dim(TAX) #1405 6 (doubletons removed)

#Create phyloseq object containing taxonomy, metadata, and otu table
phyloseq.FS9 <- phyloseq(OTU2, SAM, TAX) #prune out any OTUs that have total to 0 in all samples
phyloseq.FS9
#otu_table()   OTU Table:         [ 901 taxa and 45 samples ]
#sample_data() Sample Data:       [ 45 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 901 taxa by 6 taxonomic ranks ]

#save(phyloseq.FS9, file="phyloseq.FS9.doubleton.RData") #Use for FS9_alpha_beta_diversity.R when loading phyloseq.FS9 object

#Distance calculation
phyloseq.vegdistd7 <- vegdist(phyloseq.FS9@otu_table, method="bray")
dispersiond7 <- betadisper(phyloseq.vegdistd7, phyloseq.FS9@sam_data$All, type = c("median"))
anovabetad7 <- anova(dispersiond7)
anovabetad7
#Analysis of Variance Table

#Response: Distances
#             Df  Sum Sq    Mean Sq    F value  Pr(>F)
#Groups       5   0.12426   0.0248510  3.6715   0.008098 **
#Residuals    39  0.26398   0.0067686  


Tukeyd7 <- TukeyHSD(dispersiond7)
Tukeyd7

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = distances ~ group, data = df)

#$group
#                         diff        lwr        upr        p adj
#D11_NONINFnm-D11_INFnm    -0.045312513 -0.19456625 0.10394122 0.9417282
#D14_INFnm-D11_INFnm       -0.028293126 -0.15557714 0.09899089 0.9846537
#D14_NONINFnm-D11_INFnm    -0.024555717 -0.15767253 0.10856109 0.9934365
#D7_INFnm-D11_INFnm         0.113740456 -0.02339082 0.25087174 0.1535795
#D7_NONINFnm-D11_INFnm      0.041286895 -0.08862181 0.17119560 0.9299807
#D14_INFnm-D11_NONINFnm     0.017019387 -0.11798570 0.15202447 0.9989192
#D14_NONINFnm-D11_NONINFnm  0.020756795 -0.11976095 0.16127454 0.9976895
#D7_INFnm-D11_NONINFnm      0.159052969  0.01472646 0.30337947 0.0234448
#D7_NONINFnm-D11_NONINFnm   0.086599408 -0.05088304 0.22408186 0.4249535
#D14_NONINFnm-D14_INFnm     0.003737408 -0.11318042 0.12065524 0.9999988
#D7_INFnm-D14_INFnm         0.142033582  0.02056473 0.26350244 0.0138279
#D7_NONINFnm-D14_INFnm      0.069580021 -0.04367176 0.18283180 0.4527139
#D7_INFnm-D14_NONINFnm      0.138296173  0.01072836 0.26586399 0.0268908
#D7_NONINFnm-D14_NONINFnm   0.065842613 -0.05392729 0.18561252 0.5736010
#D7_NONINFnm-D7_INFnm      -0.072453561 -0.19667004 0.05176292 0.5101138


phyloseq.adonis <- as(sample_data(phyloseq.FS9), "data.frame")
set.seed(1)
adonis.FS9 <- adonis(phyloseq.vegdistd7~Day*Treatment, data=phyloseq.adonis, permutations=9999)
adonis.FS9

#Call:
#adonis(formula = phyloseq.vegdist ~ Day * Treatment, data = phyloseq.adonis,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                 Df    SumsOfSqs   MeanSqs   F.Model   R2        Pr(>F)    
# Day             2     1.1325      0.56626   3.07605   0.12833   0.0001 ***
#Treatment        1     0.2747      0.27474   1.49243   0.03113   0.1047    
#Day:Treatment    2     0.2386      0.11930   0.64807   0.02704   0.9519    
#Residuals        39    7.1794      0.18409             0.81351           
#Total            44    8.8253                          1.00000              
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1