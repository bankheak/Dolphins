# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# FUNCTIONS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/OneDrive/Documents/Homework/OSU/R Homework")


# SIMPLE-RATIO INDEX ------------------------------------------------------
#' @description This function creates an association matrix using the simple-ratio index (SRI). 
#' @param matr A binary matrix depicting individuals in the columns and groups in the rows
#' @return A square matrix in which each cell is an estimate of a dyadic social relationship, from 0 (never seen in the same group) to 1 (always seen in the same group)

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


# Null Permutations -------------------------------------------------------
#' @description shuffles binary matrices under different restrictions. 
#' @param mat A quantitative matrix
#' @param iter Number of random matrices to be created 
#' @param model Function to be chosen.
#' @param ... Further arguments from \code{permatswap} or \code{permatfull}
#' @return a list with \code{iter} random matrices
#' @details Totally restricted null model is called. Cell values are permuted restricting all features of the original matrix: column sums, row sums, matrix fill and total sum.
#' @references \code{citation("vegan")}

null <- function (mat, iter, ...){
  require(vegan)
  aux <- permatswap(mat, times=iter, method="quasiswap", fixedmar="both", shuffle="both", mtype="prab")
  return(aux$perm)
}

