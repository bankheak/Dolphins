
rm(list=ls()) 
gc()

# Load the parallel package
library(doParallel)
library(foreach)

# Read file in
gbi <- readRDS("../data/gbi.RData")

# Specify the number of nodes/workers in the cluster
num_nodes <- 4

# Create a cluster with the specified number of nodes/workers
cl <- makeCluster(num_nodes)

# Register the cluster to enable parallel processing
registerDoParallel(cl)

# Load Null function
null <- function (mat, iter, ...){
  library(vegan)
  aux <- permatswap(mat, times=iter, method="quasiswap", fixedmar="both", 
                    shuffle="both", mtype="prab")
  return(aux$perm)
}

# Load SRI function
SRI.func <-  function (matr) {
  if (any(is.na(matr))) {
    matr <- na.omit(matr)
    cat("The data matrix contains NA, and have been removed.\n")
  }
  matr1 = matr
  N <- nrow(matr1)
  matr1[matr1 > 1] <- 1
  n <- apply(matr1, 2, sum)
  tmatr <- t(matr1)
  df <- as.matrix(t(matr))
  a <- df %*% t(df) # Dyad in same group
  b <- df %*% (1 - t(df)) # A present, B absent
  c <- (1 - df) %*% t(df) # A absent, B present
  d <- ncol(df) - a - b - c # Double absent
  Dice <- data.frame()
  for (i in 1:nrow(a)) {
    for (j in 1:ncol(a)) {
      Dice[i, j] <- a[i, j]/(a[i, j] + b[i, j] + c[i, j])
    }
  }
  rownames(Dice)=colnames(Dice)=colnames(matr)
  Dice
}

# Parallel processing
reps <- 1000
nF <- lapply(gbi, function (df) {null(df, iter=reps)})

#' Calculate the association and CV for each of the 1000 permuted matrices to
#' create null distribution
cv_years <- lapply(nF, function(df) {
  cv_null <- foreach(i = 1:reps, .combine = c, .export = c('SRI.func')) %dopar% {
    sri_null <- as.matrix(SRI.func(df[[i]]))
    (sd(sri_null) / mean(sri_null)) * 100
  }
  cv_null
})

# Stop the cluster
stopCluster(cl)

# Write the results to a file
saveRDS(cv_years, "../code/cv_years.RData")
