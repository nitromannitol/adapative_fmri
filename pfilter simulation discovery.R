source('pfilter.R')
data <- read.table("data.txt")
P = data[, 1]
alphas = c(0.2, 0.3)
n = 10000
groups = matrix(ncol = 2, nrow = n)
groups[1:n, 1] = c(1:n)
for (j in 0:9){
  groups[(1 + j*1000):(1000 + j*1000), 2] = rep(rep((1 + j*10):(10 + j*10), each = 10), 10)
}
discoveries = pfilter(P, alphas, groups)
write.table(discoveries, "pfilter simulation discovery.txt", row.names = FALSE, col.names = FALSE)