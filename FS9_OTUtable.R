#########################################
#FS9 OTU Table
#Kathy Mou

#NOTES: 
#This code generates an OTU table

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory (either on desktop or network drive)
setwd("C:/Users/Kathy.Mou/Desktop/FSEP_Projects/FS9/FS9_RWorkspace")

#Load library packages
library(tidyverse)


#Load image file
load("FS9_OTU_table.RData")


#Save image file
save.image(file="FS9_OTU_table.RData")

#########################################
#Files needed:
#stability.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#stability.outsingletons.abund.opti_mcc.0.03.subsample.shared

#Edit taxonomy file
#Save "stability.outsingletons.abund.opti_mcc.0.03.cons.taxonomy" as a csv file in a spreadsheet editor
taxonomy <- read.csv("stability.outsingletons.abund.opti_mcc.0.03.cons.taxonomy.csv")
taxonomy
taxonomy$Taxonomy <- gsub('*\\(.*?\\) *', '', taxonomy$Taxonomy) #Remove (100) from "Taxonomy" column
taxonomy
taxonomy[1:6,3] #Show Taxonomy column rows 1 through 6
write.csv(taxonomy, file = "FS9.taxonomy.csv") #In excel, remove the "Size" and numbered rownames

#Edit subsample.shared file
#Save "stability.outsingletons.abund.opti_mcc.0.03.subsample.shared" as a csv file in a spreadsheet editor; remove "label", and "numOtus" rows
shared <- read.csv("stability.outsingletons.abund.opti_mcc.0.03.subsample.shared.csv")
head(shared)
shared[1:6,1]
shared <- t(shared) #Transpose
head(shared)
shared[1:6,1]
write.csv(shared, file = 'FS9.subsample.shared.csv') #Open in a spreadsheet editor and remove the "V*" row

#Read new subsample.shared file and edit column name
shared <- read.csv("FS9.subsample.shared.csv")
head(shared)
colnames(shared) [1] <- "OTU" #Rename first column to "OTU"

#Read new taxonomy file and merge with 'shared' dataframe
taxonomy <- read.csv("FS9.taxonomy.csv")
OTUtable <- merge(shared, taxonomy, by.x ="OTU", by.y = "OTU") #Merge "FS9.subsample.shared.csv" and "FS9.taxonomy.csv via OTU"
head(OTUtable)
nrow(OTUtable) #1583
ncol(OTUtable) #170
write.csv(OTUtable, file= "FS9.OTUtable.csv") #In a spreadsheet editor, open "FS9.OTUtable.csv" and remove first column (unnecessary)

#check OTU table
OTUtable <-read.csv("FS9.OTUtable.csv")
head(OTUtable)
