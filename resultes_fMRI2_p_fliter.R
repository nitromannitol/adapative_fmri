source("p-filter.R")
table = read.table("Desktop/308/project/fMRI_2.txt", header = TRUE)
head(table)
P<-table$pvalue
voxel<-table$voxel
ROI<-table$ROI
alphas<-c(0.05, 0.05)
group<-cbind(voxel, ROI)
Discoveries_pfilter<-pfilter(P, alphas, group)
Discoveries_pfilter
table_results<-table
table_results$discoveries<-Discoveries_pfilter   ##table added discovery of signal
write.table(table_results, "results_fMRI2_pfilter.txt")
