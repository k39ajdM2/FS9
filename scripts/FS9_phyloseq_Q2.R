#########################################
#FS9 16S phyloseq - Infected (days 4, 7): inject vs feed vs nm
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
#R version 4.0.2 (2020-06-22)

#########################################

#Read files for metadata and OTU table
meta <- read.csv("./data/FS9_metadata_INF_InjectFeedNM.csv", row.names = 1)
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
dim(otu.meta) #47 1410 (doubletons removed)
rownames(otu.meta) <- otu.meta[,1] #Set column 1 as rownames for 'otu.meta'
class(otu.meta) #dataframe
otu.meta <- otu.meta[,-1]
rownames(otu.meta)

#Added an "All" column (combines 'Treatment' and 'Day' values) in new 'otu.meta2'
otu.meta2<- cbind(otu.meta) #Make second copy of otu.meta to use to include "All" column
colnames(otu.meta2)
otu.meta2$All <- with(otu.meta2, paste0(Day, sep="_", Treatment)) #Create "All" column with Day and Treatment combined
head(otu.meta2)
dim(otu.meta2) #47 1410 (doubletons removed)
head(otu.meta2[,1405:1410])
head(otu.meta2[,1:10])
otu.meta2<- otu.meta2[,c(1:4,1410,5:1409)] #Reorder columns to have "All" column after "Treatment" column
#write.csv(otu.meta2, file="FS9.otu.meta_all.doubleton.csv")

#Pull out metadata from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
meta2 <- otu.meta2[,c(1:5)] #Take columns 1-5 ("Day" to "All") from 'otu.meta2' to make 'meta2'
head(meta2)
dim(meta2) #47 5 (doubletons removed)

#Create SAM metadata table phyloseq object
SAM = sample_data(meta2, errorIfNULL = TRUE)
head(SAM)
dim(SAM) #47 5

#Pull out OTU data from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
tail(otu.meta2[,1405:1410])
otu2 <- otu.meta2[,c(6:1410)] #Select OTU columns to create 'otu2' dataframe
head(otu2[,1:10])
dim(otu2) #47 1405 (doubletons removed)
otu2.trans <- t(otu2) #Transpose otu.all to have OTUs as rownames, sample names as column names
head(otu2.trans[,1:10])
dim(otu2.trans) #1405 47 (doubletons removed)

#Merge 'tax' back into 'otu2.trans' for correct format and taxons
head(tax)
otu2.tax <- merge(otu2.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu2.tax) #1405 50 (doubletons removed)
head(otu2.tax[,1:10])
head(otu2.tax[,45:50])
row.names(otu2.tax) <- otu2.tax[,1] #Set first row of 'otu2.tax' as rownames
head(otu2.tax[,1:5])
otu2.tax <- otu2.tax[,-1] #Remove first row aka extraneous "Row.names" column from 'otu2.tax'
head(otu2.tax[,1:5])

#Split again
dim(otu2.tax) #1405 49 (doubletons removed)
head(otu2.tax[,44:49])
otu2.notax <- otu2.tax[,1:47] #take rows 1-47 to make new dataframe 'otu2.notax' (48 is "delete" column, 49 is "Taxonomy" column)
dim(otu2.notax) #1405 47 (doubletons removed)
head(otu2.notax[,1:5])
head(otu2.notax[,40:47])
class(otu2.notax) #dataframe
otu2.notax <- as.matrix(otu2.notax) #turn 'otu2.notax' into a matrix class
class(otu2.notax) #matrix array
otu2.notax
otu2.notax.trans <- t(otu2.notax)
head(otu2.notax.trans[,1:10])

#Create OTU table phyloseq object
OTU = otu_table(otu2.notax.trans, taxa_are_rows = FALSE)
head(OTU[,1:10])
dim(OTU) #47 1405 (doubletons removed)
class(OTU)
OTU2 <- prune_taxa(taxa_sums(OTU) > 0, OTU)
class(OTU2)
head(OTU2[,1:10])
taxa_sums(OTU2)
dim(OTU2) #47 866 (doubletons removed)

#Edit taxonomy
dim(otu2.tax) #1405 49 (doubletons removed)
head(otu2.tax[,40:49])
tax2 <- separate(data = otu2.tax, 
                 col = Taxonomy, 
                 into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#"separate" function separates "Taxonomy" column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax2) #notice that Species column is blank
dim(tax2) #1405 55 (doubletons removed)
head(tax2[,49:55])
tax2.kg <- tax2[,49:54] #Keep only taxonomy columns "Kingdom" to "Genus"
head(tax2.kg)
dim(tax2.kg) #1405 6 (doubletons removed)
class(tax2.kg) #dataframe
tax2.kg <- as.matrix(tax2.kg)
class(tax2.kg) #matrix array
tax2.kg

#Create TAX taxonomy table phyloseq object
TAX = tax_table(tax2.kg)
head(TAX)
dim(TAX) #1405 6 (doubletons removed)

#Create phyloseq object containing taxonomy, metadata, and otu table
phyloseq.FS9 <- phyloseq(OTU2, SAM, TAX) #prune out any OTUs that have total to 0 in all samples
phyloseq.FS9
#otu_table()   OTU Table:         [ 866 taxa and 47 samples ]
#sample_data() Sample Data:       [ 47 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 866 taxa by 6 taxonomic ranks ]

save(phyloseq.FS9, file="phyloseq.FS9.doubleton.RData") #Use for FS9_alpha_beta_diversity.R when loading phyloseq.FS9 object

###############################################################################

#Distance calculation
phyloseq.vegdistq2 <- vegdist(phyloseq.FS9@otu_table, method="bray") 
dispersionq2 <- betadisper(phyloseq.vegdistq2, phyloseq.FS9@sam_data$All, type = c("median"))
anovabetaq2 <- anova(dispersionq2)
anovabetaq2
#Analysis of Variance Table
#Response: Distances
#           Df  Sum Sq    Mean Sq     F value Pr(>F)
#Groups     5   0.04097   0.0081941   0.7031  0.6243
#Residuals 41   0.47781   0.0116539


Tukeyq2 <- TukeyHSD(dispersionq2)
Tukeyq2

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = distances ~ group, data = df)

#$group
#                           diff         lwr        upr         p adj
#D4_INFinject-D4_INFfeed    0.0208560458 -0.18136536 0.22307746 0.9995987
#D4_INFnm-D4_INFfeed        0.0215626427 -0.18669655 0.22982184 0.9995910
#D7_INFfeed-D4_INFfeed      0.0250819916 -0.16579071 0.21595470 0.9987006
#D7_INFinject-D4_INFfeed   -0.0533979794 -0.24427068 0.13747473 0.9589152
#D7_INFnm-D4_INFfeed       -0.0067167257 -0.19758943 0.18415598 0.9999980
#D4_INFnm-D4_INFinject      0.0007065969 -0.17879026 0.18020345 1.0000000
#D7_INFfeed-D4_INFinject    0.0042259459 -0.15476971 0.16322160 0.9999995
#D7_INFinject-D4_INFinject -0.0742540251 -0.23324968 0.08474163 0.7292512
#D7_INFnm-D4_INFinject     -0.0275727714 -0.18656842 0.13142288 0.9951486
#D7_INFfeed-D4_INFnm        0.0035193490 -0.16308801 0.17012671 0.9999998
#D7_INFinject-D4_INFnm     -0.0749606220 -0.24156798 0.09164673 0.7586490
#D7_INFnm-D4_INFnm         -0.0282793683 -0.19488672 0.13832799 0.9956123
#D7_INFinject-D7_INFfeed   -0.0784799710 -0.22276617 0.06580623 0.5869497
#D7_INFnm-D7_INFfeed       -0.0317987173 -0.17608492 0.11248749 0.9854368
#D7_INFnm-D7_INFinject      0.0466812537 -0.09760495 0.19096746 0.9257214

phyloseq.adonis <- as(sample_data(phyloseq.FS9), "data.frame")
set.seed(1)
adonis.FS9 <- adonis(phyloseq.vegdist~Day*Treatment, data=phyloseq.adonis, permutations=9999)
adonis.FS9

#Call:
#adonis(formula = phyloseq.vegdist ~ Day * Treatment, data = phyloseq.adonis,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                 Df    SumsOfSqs   MeanSqs   F.Model   R2        Pr(>F)    
#Day              1     0.4396      0.43962   2.8200    0.05725   0.0023 **
#Treatment        2     0.5273      0.26365   1.6912    0.06867   0.0236 * 
#Day:Treatment    2     0.3200      0.15998   1.0262    0.04167   0.3989   
#Residuals        41    6.3916      0.15589             0.83241          
#Total            46    7.6785                          1.00000  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1