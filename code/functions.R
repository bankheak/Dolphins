# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# FUNCTIONS
###########################################################################

# load necessary packages
library(asnipe)
library(Matrix)
library(matrixStats)
library(foreach)
library(doParallel)
library(tcltk)
library(survival)
library(combinat)


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


# EDGE LIST -----------------------------------------------------------------
matrix_to_edgelist=function(mat, rawdata, idnodes){
  
  if (rawdata == TRUE){
    ID=colnames(mat)
    mat=HWI(mat)
    colnames(mat)=rownames(mat)
  } else {
    ID=colnames(mat)
    rownames(mat)=colnames(mat)=c(1:nrow(mat))
  }
  
  if (idnodes == FALSE){
    vi=vector(mode="numeric", length=nrow(mat)^2)
    vi=rep(rownames(mat), nrow(mat))
    vi=as.numeric(vi)
    
    vj=vector(mode="numeric", length=ncol(mat)^2)
    v0=colnames(mat)
    v00=matrix(0,ncol(mat), ncol(mat))
    for(i in 1:ncol(mat)){
      v00[,i]=rep(v0[i], ncol(mat))
    }
    vj=as.vector(v00)
    vj=as.numeric(vj)
    
    m0=data.matrix(mat)
    dimnames(m0)=NULL
    vw=as.vector(m0)
    
    edgelist=data.frame(vi, vj, vw)
    
    for(i in 1:(nrow(edgelist))){
      if(edgelist[i,1] == edgelist[i,2]) edgelist[i,3]=NA
    }
    edge=edgelist[!is.na(edgelist$vw),]
    
    for(i in 1:(nrow(edge))){
      if(edge[i,3] == 0.00) edge[i,3]=NA
    }
    ed=edge[!is.na(edge$vw),]
    
    eltnet=data.matrix(ed)
    as.tnet(eltnet)
    return(eltnet)
  }
  
  if (idnodes == TRUE){
    nodes = data.frame(colnames(mat), ID)
    names(nodes)=c("tnet code", "real ID")
    print(nodes)}
}


# SIMULARITY INDEX ------------------------------------------------------
sim.func<-  function (matr) {
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
  a[a > 1] <- 1
  a
}
