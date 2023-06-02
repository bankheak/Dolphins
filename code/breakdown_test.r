# Slowly build up the code to work in the cluster
require(doParallel) # registerDoParallel
require(vegan)
# Null function
null <- function (mat, iter, ...){
  aux <- permatswap(mat, times=iter, method="quasiswap", fixedmar="both", shuffle="both", mtype="prab")
  return(aux$perm)
}

# Read file in
gbi<-  read.csv("../data/gbi.csv")

#  Create 1000 random group-by-individual binary matrices
reps<- 1000
registerDoParallel(10)
nF <- null(gbi, iter=reps)

stopImplicitCluster()

write.csv(nF, "nF.csv")