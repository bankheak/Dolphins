# Load all necessary packages
# Run multiple cores for faster computing
require(doParallel) # registerDoParallel
require(foreach)

# Load function
SRI.func<-  function (matr) {
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

# Specify the number of nodes/workers in the cluster
num_nodes <- 4

# Create a cluster with the specified number of nodes/workers
cl <- makeCluster(num_nodes)

# Register the cluster to enable parallel processing
registerDoParallel(cl)

# Read file in
nF<-  readRDS("../code/nF.RData")

#' Calculate the association and CV for each of the 1000 permuted matrices to
#' create null distribution
reps <- 1000
cv_null <- rep(NA,reps)

cv_null_years <- list()
for (j in 1:22) {
  cv_null_years[[j]] <- foreach(i = 1:reps, 
                     .combine = c) %dopar% { 
                       sri_null = as.matrix(SRI.func(nF[[i]]))
                       cv_null[i] <- ( sd(sri_null) / mean(sri_null) ) * 100} 
}

# Stop the cluster
stopCluster(cl)

# Write the results to a file
write.csv(cv_null, "cv_null.csv", row.names = FALSE)
