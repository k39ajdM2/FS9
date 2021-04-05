#########################################
#FS9 16S phyloseq - Infected (days 7, 11, 14): inject vs feed vs nm
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
dim(otu.meta) #72 1410 (doubletons removed)
rownames(otu.meta) <- otu.meta[,1] #Set column 1 as rownames for 'otu.meta'
class(otu.meta) #dataframe
otu.meta <- otu.meta[,-1]
rownames(otu.meta)

#Added an "All" column (combines 'Treatment' and 'Day' values) in new 'otu.meta2'
otu.meta2<- cbind(otu.meta) #Make second copy of otu.meta to use to include "All" column
colnames(otu.meta2)
otu.meta2$All <- with(otu.meta2, paste0(Day, sep="_", Treatment)) #Create "All" column with Day and Treatment combined
head(otu.meta2)
dim(otu.meta2) #72 1410 (doubletons removed)
head(otu.meta2[,1405:1410])
head(otu.meta2[,1:10])
otu.meta2<- otu.meta2[,c(1:4,1410,5:1409)] #Reorder columns to have "All" column after "Treatment" column
#write.csv(otu.meta2, file="FS9.otu.meta_all.doubleton.csv")

#Pull out metadata from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
meta2 <- otu.meta2[,c(1:5)] #Take columns 1-5 ("Day" to "All") from 'otu.meta2' to make 'meta2'
head(meta2)
dim(meta2) #72 5 (doubletons removed)

#Create SAM metadata table phyloseq object
SAM = sample_data(meta2, errorIfNULL = TRUE)
head(SAM)
dim(SAM) #72 5

#Pull out OTU data from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
tail(otu.meta2[,1405:1410])
otu2 <- otu.meta2[,c(6:1410)] #Select OTU columns to create 'otu2' dataframe
head(otu2[,1:10])
dim(otu2) #72 1405 (doubletons removed)
otu2.trans <- t(otu2) #Transpose otu.all to have OTUs as rownames, sample names as column names
head(otu2.trans[,1:10])
dim(otu2.trans) #1405 72 (doubletons removed)

#Merge 'tax' back into 'otu2.trans' for correct format and taxons
head(tax)
otu2.tax <- merge(otu2.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu2.tax) #1405 75 (doubletons removed)
head(otu2.tax[,1:10])
head(otu2.tax[,70:75])
row.names(otu2.tax) <- otu2.tax[,1] #Set first row of 'otu2.tax' as rownames
head(otu2.tax[,1:5])
otu2.tax <- otu2.tax[,-1] #Remove first row aka extraneous "Row.names" column from 'otu2.tax'
head(otu2.tax[,1:5])

#Split again
dim(otu2.tax) #1405 74 (doubletons removed)
head(otu2.tax[,70:74])
otu2.notax <- otu2.tax[,1:72] #take rows 1-72 to make new dataframe 'otu2.notax' (73 is "delete" column, 74 is "Taxonomy" column)
dim(otu2.notax) #1405 72 (doubletons removed)
head(otu2.notax[,1:5])
head(otu2.notax[,70:72])
class(otu2.notax) #dataframe
otu2.notax <- as.matrix(otu2.notax) #turn 'otu2.notax' into a matrix class
class(otu2.notax) #matrix array
otu2.notax
otu2.notax.trans <- t(otu2.notax)
head(otu2.notax.trans[,1:10])

#Create OTU table phyloseq object
OTU = otu_table(otu2.notax.trans, taxa_are_rows = FALSE)
head(OTU[,1:10])
dim(OTU) #72 1405 (doubletons removed)
class(OTU)
OTU2 <- prune_taxa(taxa_sums(OTU) > 0, OTU)
class(OTU2)
head(OTU2[,1:10])
taxa_sums(OTU2)
dim(OTU2) #72 1042 (doubletons removed)

#Edit taxonomy
dim(otu2.tax) #1405 72 (doubletons removed)
head(otu2.tax[,70:72])
tax2 <- separate(data = otu2.tax, 
                 col = Taxonomy, 
                 into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#"separate" function separates "Taxonomy" column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax2) #notice that Species column is blank
dim(tax2) #1405 80 (doubletons removed)
head(tax2[,70:80])
tax2.kg <- tax2[,74:79] #Keep only taxonomy columns "Kingdom" to "Genus"
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
#otu_table()   OTU Table:         [ 1042 taxa and 72 samples ]
#sample_data() Sample Data:       [ 72 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 1042 taxa by 6 taxonomic ranks ]

#save(phyloseq.FS9, file="phyloseq.FS9.doubleton.RData") #Use for FS9_alpha_beta_diversity.R when loading phyloseq.FS9 object

#Distance calculation
phyloseq.vegdistq2 <- vegdist(phyloseq.FS9@otu_table, method="bray") 
dispersionq2 <- betadisper(phyloseq.vegdistq2, phyloseq.FS9@sam_data$All, type = c("median"))
anovabetaq2 <- anova(dispersionq2)
anovabetaq2
#Analysis of Variance Table
#Response: Distances
#           Df  Sum Sq    Mean Sq     F value Pr(>F)
#Groups     8   0.18789   0.0234858   2.372   0.02666 *
#Residuals 63   0.62377   0.0099012 


Tukeyq2 <- TukeyHSD(dispersionq2)
Tukeyq2

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = distances ~ group, data = df)

#$group
#                           diff         lwr        upr         p adj
#D11_INFinject-D11_INFfeed    0.021097848 -0.179228846 0.22142454 0.9999943
#D11_INFnm-D11_INFfeed        0.022275376 -0.184032533 0.22858329 0.9999931
#D14_INFfeed-D11_INFfeed      0.025840283 -0.163244039 0.21492460 0.9999572
#D14_INFinject-D11_INFfeed   -0.052660806 -0.241745128 0.13642352 0.9925093
#D14_INFnm-D11_INFfeed       -0.006161034 -0.195245356 0.18292329 1.0000000
#D7_INFfeed-D11_INFfeed       0.018948922 -0.195452546 0.23335039 0.9999985
#D7_INFinject-D11_INFfeed     0.072232306 -0.110512183 0.25497679 0.9363892
#D7_INFnm-D11_INFfeed         0.136023758 -0.064302937 0.33635045 0.4303410
#D11_INFnm-D11_INFinject      0.001177528 -0.176637528 0.17899258 1.0000000
#D14_INFfeed-D11_INFinject    0.004742435 -0.152763506 0.16224838 1.0000000
#D14_INFinject-D11_INFinject -0.073758655 -0.231264595 0.08374729 0.8494279
#D14_INFnm-D11_INFinject     -0.027258882 -0.184764823 0.13024706 0.9997448
#D7_INFfeed-D11_INFinject    -0.002148926 -0.189293870 0.18499602 1.0000000
#D7_INFinject-D11_INFinject   0.051134457 -0.098701415 0.20097033 0.9727292
#D7_INFnm-D11_INFinject       0.114925910 -0.055913270 0.28576509 0.4432102
#D14_INFfeed-D11_INFnm        0.003564907 -0.161481420 0.16861123 1.0000000
#D14_INFinject-D11_INFnm     -0.074936182 -0.239982509 0.09011015 0.8699774
#D14_INFnm-D11_INFnm         -0.028436410 -0.193482737 0.13660992 0.9997531
#D7_INFfeed-D11_INFnm        -0.003326454 -0.196860427 0.19020752 1.0000000
#D7_INFinject-D11_INFnm       0.049956930 -0.107786398 0.20770026 0.9827960
#D7_INFnm-D11_INFnm           0.113748382 -0.064066673 0.29156344 0.5125723
#D14_INFinject-D14_INFfeed   -0.078501089 -0.221435401 0.06443322 0.7048042
#D14_INFnm-D14_INFfeed       -0.032001317 -0.174935629 0.11093300 0.9983334
#D7_INFfeed-D14_INFfeed      -0.006891361 -0.181949426 0.16816671 1.0000000
#D7_INFinject-D14_INFfeed     0.046392023 -0.088043407 0.18082745 0.9708399
#D7_INFnm-D14_INFfeed         0.110183475 -0.047322466 0.26768942 0.3894045
#D14_INFnm-D14_INFinject      0.046499772 -0.096434540 0.18943408 0.9796540
#D7_INFfeed-D14_INFinject     0.071609729 -0.103448337 0.24666779 0.9234443
#D7_INFinject-D14_INFinject   0.124893112 -0.009542318 0.25932854 0.0887652
#D7_INFnm-D14_INFinject       0.188684564  0.031178623 0.34619051 0.0080540
#D7_INFfeed-D14_INFnm         0.025109956 -0.149948110 0.20016802 0.9999380
#D7_INFinject-D14_INFnm       0.078393339 -0.056042090 0.21282877 0.6342041
#D7_INFnm-D14_INFnm           0.142184792 -0.015321149 0.29969073 0.1084185
#D7_INFinject-D7_INFfeed      0.053283383 -0.114906966 0.22147373 0.9827596
#D7_INFnm-D7_INFfeed          0.117074836 -0.070070109 0.30421978 0.5428794
#D7_INFnm-D7_INFinject        0.063791452 -0.086044420 0.21362732 0.9056417

phyloseq.adonis <- as(sample_data(phyloseq.FS9), "data.frame")
set.seed(1)
adonis.FS9 <- adonis(phyloseq.vegdistq2~Day*Treatment, data=phyloseq.adonis, permutations=9999)
adonis.FS9

#Call:
#adonis(formula = phyloseq.vegdist ~ Day * Treatment, data = phyloseq.adonis,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                 Df    SumsOfSqs   MeanSqs   F.Model   R2        Pr(>F)    
#Day              2     1.4053      0.70263   3.9052    0.10012   0.0001 ***
#Treatment        2     0.6112      0.30561   1.6986    0.04355   0.0271 *  
#Day:Treatment    4     0.6842      0.17105   0.9507    0.04875   0.5448    
#Residuals        63    11.3351     0.17992             0.80759           
#Total            71    14.0358                         1.00000  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1