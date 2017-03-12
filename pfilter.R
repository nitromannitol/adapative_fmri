pfilter = function(P,alphas,groups){
	# P in [0,1]^n = vector of p-values
	# alphas in [0,1]^M = vector of target FDR levels
	# groups is a n-by-M matrix; 
	#	groups[i,m] = which group does P[i] belong to,
	#		for the m-th grouping
	
	n = length(P)
	M = length(alphas)
	G = apply(groups,2,max) # G[m] = # groups, for grouping m
	
	
	Simes = list()
	for(m in 1:M){
		Simes[[m]]=rep(0,G[m])
		for(g in 1:G[m]){
			group = which(groups[,m]==g)
			Simes[[m]][g]=min(sort(P[group])*length(group)/(1:length(group)))
		}
	}
	
		
	# initialize
	thresh = alphas
	Sh = 1:n
	for(m in 1:M){
		pass_Simes_m = which(is.element(groups[,m],which(Simes[[m]]<=thresh[m])))
		Sh = intersect(Sh,pass_Simes_m)
	}
	done = FALSE


	while(!done){
		thresh_old = thresh
		for(m in 1:M){
			# which groups, for the m-th grouping, 
			#	have any potential discoveries?
			Shm = sort(unique(groups[Sh,m]))

			# run BH on Simes[[m]], constraining to Shm
			Pvals_m = rep(1.01,G[m]); # >1 for groups not in Dm
			Pvals_m[Shm] = Simes[[m]][Shm]
			khatm = max(0,which(sort(Pvals_m)<=(1:G[m])/G[m]*alphas[m]))
			thresh[m] = khatm/G[m]*alphas[m]
			Sh = intersect(Sh,
				which(is.element(groups[,m],which(Pvals_m<=thresh[m]))))
		}
		if(all(thresh_old==thresh)){done = TRUE}
	}
	
	Sh_temp = Sh;
	Sh = rep(0,n); Sh[Sh_temp] = 1
	Sh
	
}