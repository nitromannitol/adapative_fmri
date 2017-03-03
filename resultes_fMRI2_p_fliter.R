source("p-filter.R")
table = read.table("Desktop/308/project/fMRI_2.txt", header = TRUE)
head(table)
P<-table$pvalue
voxel<-table$voxel
ROI<-table$ROI
alphas<-c(0.05, 0.05)
group<-cbind(voxel, ROI)
Discoveries_pfilter<-pfilter(P, alphas, group)
Discoveries_pfiler
table_results<-cbind(table[, 1:4], Discoveries_pfilter) 
head(table_results)
