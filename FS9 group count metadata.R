library(dplyr)

FS9 <- read.csv("FS9 final group count.csv")
FS9$All <- with(FS9, paste0(FS9$Day, sep="_", FS9$Group))
distinct(FS9$All)
table(FS9$All)
write.csv(FS9, file="FS9_Final_Group_Count_with_All.csv", row.names = TRUE)
