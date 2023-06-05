# Slowly build up the code to work in the cluster
# Load all necessary packages
require(asnipe) # get_group_by_individual--Damien Farine
# Could do permutations
require(assocInd)
# Run multiple cores for faster computing
require(doParallel) # registerDoParallel
require(microbenchmark)
require(parallel)
require(foreach)
# Null function
null <- function (mat, iter, ...){
  require(vegan)
  aux <- permatswap(mat, times=iter, method="quasiswap", fixedmar="both", shuffle="both", mtype="prab")
  return(aux$perm)
}

# Read file in
gbi<-  readRDS("../data/gbi.RData")

#  Create 1000 random group-by-individual binary matrices
gbi <- gbi[[1]]
reps<- 1000
registerDoParallel(2)
nF <- null(gbi, iter=reps)

#' Calculate the association and CV for each of the 1000 permuted matrices to
#' create null distribution
cv_null <- rep(NA,reps)


    foreach(i = 1:reps, 
            .combine = c) %dopar% { 
            sri_null = as.matrix(SRI.func(nF[[i]]))
            cv_null[i] <- ( sd(sri_null) / mean(sri_null) ) * 100}

stopImplicitCluster()

# remove NAs, if any
cv_null = cv_null[!is.na(cv_null)]

saveRDS(cv_null, file = "/home/ib/bankheak/R_libs/Dolphins/code/cv_null.rds")


# Perform the computation on each chunk in parallel
nF <- foreach(chunk = array_chunks, .combine = c) %dopar% {
  #  Create 1000 random group-by-individual binary matrices
  null(gbi, iter=1000)
}