#########################################
#FS9 OTU Table
#Kathy Mou

#NOTES: 
#This code generates an OTU table

#Clear workspace and load necessary packages
rm(list=ls())

#Load library packages
library(tidyverse)

#########################################
#Files needed:
#stability.outdoubletons.abund.opti_mcc.0.03.cons.taxonomy
#stability.outdoubletons.abund.opti_mcc.0.03.subsample.shared

#Edit taxonomy file
#Save "stability.outdoubletons.abund.opti_mcc.0.03.cons.taxonomy" as a csv file in a spreadsheet editor
taxonomy <- read.csv("./data/stability.outdoubletons.abund.opti_mcc.0.03.cons.taxonomy.csv")
taxonomy
taxonomy$Taxonomy <- gsub('*\\(.*?\\) *', '', taxonomy$Taxonomy) #Remove (100) from "Taxonomy" column
taxonomy
taxonomy[1:6,3] #Show Taxonomy column rows 1 through 6
write.csv(taxonomy, file = "FS9.taxonomy.doubleton.csv") #In excel, remove the "Size" and numbered rownames

#Edit subsample.shared file
#Save "stability.outdoubletons.abund.opti_mcc.0.03.subsample.shared" as a csv file in a spreadsheet editor; remove "label", and "numOtus" rows
shared <- read.csv("./data/stability.outdoubletons.abund.opti_mcc.0.03.subsample.shared.csv")
head(shared)
shared[1:6,1]
shared <- t(shared) #Transpose
head(shared)
shared[1:6,1]
write.csv(shared, file = 'FS9.subsample.shared.doubleton.csv') #Open in a spreadsheet editor and remove the "V*" row

#Read new subsample.shared file and edit column name
shared <- read.csv("./data/FS9.subsample.shared.doubleton.csv")
head(shared)
colnames(shared) [1] <- "OTU" #Rename first column to "OTU"

#Read new taxonomy file and merge with 'shared' dataframe
taxonomy <- read.csv("./data/FS9.taxonomy.doubleton.csv")
OTUtable <- merge(shared, taxonomy, by.x ="OTU", by.y = "OTU") #Merge 'shared' and 'taxonomy' via "OTU"
head(OTUtable)
nrow(OTUtable) #1583 (singletons removed), 1405 (doubletons removed)
ncol(OTUtable) #170 (singletons removed), 169 (doubletons removed)
write.csv(OTUtable, file= "FS9.OTUtable.doubleton.csv") #In a spreadsheet editor, open "FS9.OTUtable.doubleton.csv" and remove first column (unnecessary)

#check OTU table
OTUtable <-read.csv("./data/FS9.OTUtable.csv")
head(OTUtable)
