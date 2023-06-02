# Load the parallel package
library(parallel)
library(doParallel)
library(pbapply)

# Read file in
gbi<-  read.csv("../data/gbi.csv")

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

# Create an example array
my_array <- 1:1000

# Split the array into chunks for parallel processing
array_chunks <- split(my_array, ceiling(seq_along(my_array) / num_nodes))

# Perform the computation on each chunk in parallel
nF <- foreach(chunk = array_chunks, .combine = c) %dopar% {
  #  Create 1000 random group-by-individual binary matrices
  null(gbi, iter=1000)
  
}

# Stop the cluster
stopCluster(cl)

# Print the results
print(nF)
