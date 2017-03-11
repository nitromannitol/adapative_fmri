install.packages("R.matlab")
paired_t_test = function(X1, X2){
  n = dim(X1)[1] # X1 & X2 are vectors of length n
  Diff = X1 - X2
  avg_D = colMeans(Diff)
  var_D = colSums((Diff - matrix(1, n, 1)%*%avg_D)^2) / (n-1)
  t_stat = abs(avg_D)/sqrt(var_D/n)
  p_val = (1-pt(t_stat, df = n-1))*2    
  p_val
}

fMRI_get_data_and_pvals = function(data_downloaded = FALSE){
  
  if(!require(R.matlab)){
    install.packages("R.matlab")
  }
  library("R.matlab")
  
  ### Download CMU StarPlus fMRI data (subject 04847) from web
  if (!data_downloaded){
    download.file(url = "http://www.cs.cmu.edu/afs/cs.cmu.edu/project/theo-81/www/data-starplus-04847-v7.mat", destfile = "./fMRI_data.mat")
  }
  
  ### Read "meta" and "data" from data-starplus-04847-v7.mat:
  Subject = readMat("./fMRI_data.mat")
  n = Subject$meta[[3]] # number of voxels = 4698
  voxel_coords = Subject$meta[[1]] # 3D coordinates of each voxel = n-by-3 matrix
  
  ### organize labels of each voxel into the Regions Of Interest (ROIs)
  ROI_names = rep("0", 24)  
  for (i in 1:24){
    ind = c(1:3,5:25)[i] # ROI #4 , called 'LIPG', is not assigned to any voxel
    ROI_names[i] = Subject$meta[[16]][[1+3*(ind-1)]]
  }
  
  ROI_voxels = list() # ROI_voxels[[i]] contains voxel numbers (1 - 4698) of each voxel in ROI #i
  ROI_labels = rep(0,n) # ROI_labels[i] is the ROI label of voxel #i
  for(i in 1:n){
    which_ROI = which(ROI_names == toString(Subject$meta[[17]][[i]][[1]]))
    if(length(which_ROI)>0){
      ROI_labels[i] = which_ROI
    }
  }
  
  # remove voxels which are not assigned to an ROI
  remove_voxels = which(ROI_labels==0)
  n = n - length(remove_voxels)
  ROI_labels = ROI_labels[-remove_voxels]
  voxel_coords = voxel_coords[-remove_voxels,]
  
  ### Process the activity data
  # Find which trials have Picture (P) phase first, Sentence (S) phase second
  trials = which((Subject$info[1,,] > 1) & (Subject$info[14,,] == "P"))
  
  # Get the average activity levels recorded in each of the selected trials, for the P & S phases
  P_act = S_act = matrix(0, length(trials), n)
  for (i in 1:length(trials)){
    P_act[i,] = colMeans(Subject$data[[trials[i]]][[1]][1:16,-remove_voxels]) # time 1-16 = P phase
    S_act[i,] = colMeans(Subject$data[[trials[i]]][[1]][17:32,-remove_voxels]) # time 17-32 = S phase
  }
  ### Compute p-values via a paired t-test
  pvals = matrix(paired_t_test(P_act, S_act), n, 1)  
  
  output = list()
  output$P_act = P_act
  output$S_act = S_act
  output$pvals = pvals
  output$ROI_names = ROI_names
  output$ROI_labels = ROI_labels
  output$voxel_coords = voxel_coords
  output
}
output = fMRI_get_data_and_pvals()
p_values<-output$pvals
voxel_cor<-output$voxel_coords
roi<-output$ROI_labels
fmri<-cbind(1:4691, voxel_cor, roi, p_values)
write(t(fmri), "fMRI_2.txt", ncolumns = 6 )
 
