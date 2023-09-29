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
library(rgeos)
library(raster)


# SIMPLE-RATIO INDEX ------------------------------------------------------
#' @description This function creates an association matrix using the simple-ratio index (SRI). 
#' @param matr A binary matrix depicting individuals in the columns and groups in the rows
#' @return A square matrix in which each cell is an estimate of a dyadic social relationship, from 0 (never seen in the same group) to 1 (always seen in the same group)

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


# NULL PERMUTATIONS -------------------------------------------------------
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

# HRO ------------------------------------------------------
#' @description This function calculateS home range overlap using the formula: 
#' HRO = (Rij/Ri) * (Rij/Rj), where Rij is the overlap in range between individuals i and j, 
#' Ri is the range of individual i, and Rj is the range of individual j. 
#' @param kernel.poly home range polygons for each individual in a SpatialPolygonsDataFrame
#' @return A square matrix in which each cell is an estimate of dyadic homerange overlap

HRO <- function(dolph.kernel.poly) {
  
  # Empty matrix to store HRO values for each pair of individuals
  n_individuals <- length(dolph.kernel.poly@data$id)
  HRO_matrix <- matrix(NA, nrow = n_individuals, ncol = n_individuals, 
                       dimnames = list(dolph.kernel.poly@data$id, dolph.kernel.poly@data$id))
  
  # Loop through each pair of individuals to calculate HRO
  for (i in 1:n_individuals) {
    for (j in (i+1):n_individuals) {
      # Get the home range polygons for individuals i and j
      poly_i <- dolph.kernel.poly@polygons[[i]]
      poly_j <- dolph.kernel.poly@polygons[[j]]
      
      # Calculate the intersection area between the home ranges (Rij)
      intersection_area <- function(poly1, poly2) {
        if (is.null(gIntersects(poly1, poly2))) {
          return(0)  # No intersection, return 0
        } else {
          intersect_poly <- gIntersection(poly1, poly2)
          return(gArea(intersect_poly))
        }
      }
      
      intersect_area <- intersection_area(poly_i, poly_j)
      
      # Calculate the home range area for individuals i and j (Ri and Rj)
      home_range_area_i <- st_area(poly_i)
      home_range_area_j <- st_area(poly_j)
      
      # Calculate HRO using the formula: HRO = (Rij/Ri) * (Rij/Rj)
      HRO_value <- (intersect_area / home_range_area_i) * (intersect_area / home_range_area_j)
      
      # Store the HRO value in the matrix
      HRO_matrix[dolph.kernel.poly@data$id[i], dolph.kernel.poly@data$id[j]] <- HRO_value
      HRO_matrix[dolph.kernel.poly@data$id[j], dolph.kernel.poly@data$id[i]] <- HRO_value
    }
  }
  
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

# POPULATION TURNOVER ------------------------------------------------------
#' @description Calculates turnover of a given population across different periods of time and compare it to the null expectation using a null model approach
#' @param data Binary matrix \code{M} of periods of time in the rows vs. individuals in the columns, where \code{m_ij = 1} represent the presence of individual \code{j} in the period of time \code{i} and \code{m_ij = 0} otherwise
#' @param iter Integer, number of iterations for the randomization process
#' @param subseq Boolean. TRUE average only the beta diversity index (turnover) between subsequent periods (1-2, 2-3, 3-4 etc). FALSE averages all matrix elements (turnover among all periods).
#' @param plot Boolean. Plot histogram, empirical mean and 95%CI of the ramdom distribution. 
#' @details The proxy measure for turnover is the averaged Whittaker's beta diversity index for all the periods of time. The null model randomizes individuals into periods of time (columns), constraining their empirical capture frequency (row sums) and total number of sightings (matrix fill).
#' @value Returns the empirical turnover (average Whittaker index), the standard deviation (SD) and the 95% confidence intervals (two-tailed test). Significant results are shown by empirical values falling outside than the 97.5% CI or 2.5%. The funciton also returns the null distribution of Whittaker indices, where the red line represents the empirical value and the blue ones represent the 95% confidence intervals.
#' @author Mauricio Cantor (m.cantor@ymail.com)
#' @references Cantor et al. 2012 Disentangling social networks from spatiotemporal dynamics: the temporal structure of a dolphin society. Animal Behaviour 84 (2012) 641-651
#' @return

turnover_w <- function(data, iter=1000, subseq=FALSE, plot=FALSE){
  
  # create internal objects
  rand = data
  turno = numeric()
  result = matrix(NA, 1, 4);  colnames(result) = c("Empirical", "SD", "2.5%CI", "97.5%CI");  rownames(result) = c("Turnover")
  
  # calculate the turnover for the empirical data
  obs = betadiver(data, method="w", order=F, help=F)
  
  if(subseq==TRUE){
    obs = as.matrix(obs)
    aux = numeric()
    for(i in 2:nrow(obs)){ aux[i-1] = obs[i, i-1] }
    obs = aux
  }
  
  # randomize original data, calculate and save mean turnover for each iteration
  for(i in 1:iter){
    rand = apply(data, 2, sample)
    
    if(subseq==TRUE){
      aux2 = as.matrix(betadiver(rand, method="w", order=F, help=F))
      aux3 = numeric()
      for(j in 2:nrow(aux2)){ aux3[j-1] = aux2[j, j-1] }
      turno[i] = mean(aux3)
    } else {
      turno[i] = mean(betadiver(rand, method="w", order=F, help=F))
    }
  }
  
  # print result
  result[,1:2] = c(mean(obs), sd(obs))
  result[,3:4] = quantile(turno, probs=c(0.025,0.975), type=2)
  
  # plot
  if(plot==TRUE){
    hist(turno, xlab="Random turnover", main=NULL, xlim=c(result[,3]-0.03, result[,4]+0.03))
    abline(v=mean(obs), col="red")         # empirical
    abline(v=mean(result[,4]), col="blue") # 2.5% CI
    abline(v=mean(result[,3]), col="blue") # 97.5% CI
  }
  
  return(result)
}



