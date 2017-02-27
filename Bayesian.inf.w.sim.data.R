# Artificial data sets
t<-5; n<-100
loc <- rep(1:n, t)
time <- rep(1:t,each=n)
p.values.all.null <- rbeta(t*n,1,1)
data1 <- cbind(location, time, p.values.all.null)
colnames(data1) <- c("loc", "time", "p")

alphas = seq(0.1,0.5,0.1)
data.hhs <- list()
for(i in 1:length(alphas)){
  p.values.half.null <- rbeta(n/2,1,1)
  p.values.half.signal <- rbeta(n/2,alphas[i],1)
  p <- rep(c(p.values.half.null, p.values.half.signal), t) + abs(rnorm(t*n,0,0.01))
  data.hhs[[i]]<- cbind(location, time, p)
}
#TODO: can simulate more artificial data sets
#relax the condition of right parameter one from beta distributions
#generate p values from non-beta distributions
#etc



# METHOD 1:
# Assumes beta distribution, right one parameter is always one
# For each entry, use set of fixed radius to find spatially and temporally "close" points and 
# calculate mle of alpha.
# Use prior distribution to either select alpha to maximize P(alpha|Phi) or calculate weighted
# average of alpha

# use Euclidean distance. loc.scale is used to scale location
distance <- function(loc1, loc2, t1, t2, loc.scale){
  dist = sqrt((loc1/loc.scale - loc2/loc.scale)^2 + (t1 - t2)^2)
  return(dist)
}

# maximum likelihood estimation of alpha (given right parameter beta = 1)
alpha.MLE.right.one <- function(X){
  alpha = -length(X)/sum(log(X))
  return(alpha)
}

# density of Phi. TODO: change from dummy to something more reasonable.
dphi.dummy <- function(alpha){
  f = 0
  if(alpha > 0 & alpha<=1){
    f = 1
  } 
  return(f)
}

bayesian.inf <- function(data){
  alphas.bay = rep(0, nrow(data))
  loc.scale = ceiling((max(data[,"loc"]) - min(data[,"loc"])) / (max(data[,"time"]) - min(data[,"time"])))
  for(i in 1:nrow(data)){
    current = data[i,]
    current.distances = rep(0, nrow(data))
    for(j in 1:nrow(data)){
      temp = data[j,]
      current.distances[j] = distance(current[1],temp[1],current[2],temp[2],loc.scale)
    }
    l = max(current.distances)/10
    u = max(current.distances)
    radius = seq(l, u, (u-l)/10)
    alphas.mle = rep(0, length(radius))
    for(k in 1:length(radius)){
      selected.index = which(current.distances < radius[k])
      X = data[selected.index, "p"]
      alphas.mle[k] = alpha.MLE.right.one(X)
    }
    current.alpha = weighted.mean(alphas.mle, sapply(alphas.mle, FUN=dphi.dummy))
    alphas.bay[i] = current.alpha
  }
  return(alphas.bay)
}

bayesian.inf(data1)
bayesian.inf(data.hhs[[1]])
bayesian.inf(data.hhs[[2]])
bayesian.inf(data.hhs[[3]])
bayesian.inf(data.hhs[[4]])
bayesian.inf(data.hhs[[5]])



# METHOD 2: Spatial Dirichlet Process
# Compare to indepdendent Dirichlet Process, alpha\weight+theta from stick-breaking process all depend on S, 
# which is a space measure.
# TODO





