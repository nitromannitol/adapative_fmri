library(scatterplot3d)
data = read.table("fMRI_full.txt");
locations = data[,1:3];
ROI = data[,4];
pvalues = data[,5];
n = dim(locations)[1];


##IMPORTANT NULLS ARE 1 IN THE FILE
predictions = 1 - read.table("BDC_predictions.txt");
labels = read.table("BDC_labels.txt");
labels = labels[,1];
signal_colors = rep("blue",n)
for(i in 1:n){
  if(predictions[i,1] == 1){
    signal_colors[i] = "red"
  }
}
jpeg('/test/BDC_predictions.jpg')
scatterplot3d(locations[,1], locations[,2], locations[,3], signal_colors, xlab = "x", ylab = "y", zlab = "z");
title("Birth death chain")
dev.off()


source("pfilter.R")
alphas = 0.2; 
n = length(pvalues)
trivial = cbind(1:n,1:n)
Discoveries_BH = pfilter(pvalues,alphas,trivial) # BH procedure
signal_colors_BH = rep("blue",n)
for(i in 1:n){
  if(Discoveries_BH[i] == 1){
    signal_colors_BH[i] = "red"
  }
}
jpeg('BH_predictions.jpg')
scatterplot3d(locations[,1], locations[,2], locations[,3], signal_colors_BH, xlab = "x", ylab = "y", zlab = "z");
title("BH at alpha = 0.1")
dev.off()


alphas = c(0.2,0.3); 
n = length(pvalues)
#groups = cbind(1:n, labels)
groups = cbind(1:n,ROI)
Discoveries_pfilter= pfilter(pvalues,alphas,groups) 
signal_colors_pfilter = rep("blue",n)
for(i in 1:n){
  if(Discoveries_pfilter[i] == 1){
    signal_colors_pfilter[i] = "red"
  }
}
jpeg('pfilter_predictions2.jpg')
scatterplot3d(locations[,1], locations[,2], locations[,3], signal_colors_pfilter, xlab = "x", ylab = "y", zlab = "z");
title("pfilter with two groups at levels alpha = 0.2 and 0.3")
dev.off()
