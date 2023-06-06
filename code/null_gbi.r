# Load the parallel package
library(parallel)
library(doParallel)
library(foreach)

# Read file in
gbi<-  readRDS("../data/gbi.RData")

# Specify the number of nodes/workers in the cluster
num_nodes <- 4

# Create a cluster with the specified number of nodes/workers
cl <- makeCluster(num_nodes)

# Register the cluster to enable parallel processing
registerDoParallel(cl)

# Define a function that performs the computation on a single element
# Null function
null <- function (mat, iter, ...){
  library(vegan)
  aux <- permatswap(mat, times=iter, method="quasiswap", fixedmar="both", 
                    shuffle="both", mtype="prab")
  return(aux$perm)
}

# Parallel processing
reps <- 1000
nF <- list()
for (i in 1:22) {
  nF[[i]] <- null(gbi[[i]], iter=reps)
}

# Stop the cluster
stopCluster(cl)
stopImplicitCluster()

# Print the results
saveRDS(nF, "../code/nF.RData")
