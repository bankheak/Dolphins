# Load the parallel package
library(parallel)
library(doParallel)
library(foreach)

# Read file in
gbi<-  read.csv("../data/test_gbi.csv")

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
nF <- null(gbi, iter=1000)

# Stop the cluster
stopCluster(cl)
stopImplicitCluster()

# Print the results
write.csv(nF, "../code/nF.csv")
