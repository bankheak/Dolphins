# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# FUNCTIONS
###########################################################################

# load necessary packages
library(asnipe)
library(Matrix)
library(survival)
library(combinat)
library(matrixStats)
library(foreach)
library(doParallel)
library(tcltk)


# TOGETHER_APART----- ------------------------------------------------------
#' @description This function extracts those numbers from the sightings record and 
#' stores the two matrices in an array. These are then used as a basis to simulate 
#' an association matrix for the simulation part, representing a social network with noise. 
#' @param sightings group by individual matrix with binary entries
#' @return The output is an array with two matrices, the first one containing the number of 
#' times each dyad has been seen together and the second one containing the number of times 
#' members of a dyad have been seen apart. This can be used for the simulation on the sensitivity of NBDA. 
#' @references \code{citation(Wild and Hoppitt 2019)}

together_apart <- function(sightings){
  
  IDs <- colnames(sightings)
  num_sightings <- as.data.frame(colSums(sightings))
  
  # create a matrix that contains a measure of how many times each dyad has been seen together
  together_matrix <- matrix(-99, nrow=length(IDs), ncol=length(IDs))
  row.names(together_matrix) <- IDs # assign IDs as rownames
  colnames(together_matrix) <- IDs # assign IDs as colnames
  for (j in 1:length(IDs)){ # loop through each individual
    for (i in 1:length(IDs)){
      sub <- subset(sightings, subset=sightings[,i]==1 &sightings[,j]==1) # subset the data to where individuals i and j have been seen together
      together <- length(sub[,1]) # extract how many times this was the case
      together_matrix[j,i] <- together # save in matrix 
      together_matrix[i,j] <- together # save in matrix
    } # close inner loop
  } # close outer loop
  
  # create a matrix that contains a measure of how many times each dyad has been seen independenly of each other
  apart_matrix <- matrix(-99, nrow=length(IDs), ncol=length(IDs))
  row.names(apart_matrix) <- IDs # assign IDs as rownames
  colnames(apart_matrix) <- IDs # assign IDs as colnames
  
  for (j in 1:length(IDs)){ # loop through each individual
    for (i in 1:length(IDs)){
      sub <- subset(sightings, subset=sightings[,i]==1 &sightings[,j]!=1) # subset the data to where individual i has been seen without j
      apart_i <- length(sub[,1]) # extract how many times this was the case
      
      sub2 <- subset(sightings, subset=sightings[,i]!=1 &sightings[,j]==1) # subset the data to where individuals j has been seen without i
      apart_j <- length(sub2[,1])
      
      apart_matrix[j,i] <- apart_i+apart_j  # save in matrix 
      apart_matrix[i,j] <- apart_i+apart_j  # save in matrix
    } # close inner loop
  } # close outer loop
  
  # save the two matrices in an array with two dimensions
  together_apart_array <-  array(c(together_matrix,apart_matrix), 
                                 dim=c(length(together_matrix[,1]), length(apart_matrix[,1]), 2),
                                 dimnames=c(list(IDs),list(IDs))) # assign IDs as column and rownames
  
  
  return(together_apart_array)
  
} #close function loop


# SENSITIVITY NBDA ERROR-----------------------------------------------------------
#' @description This function simulates a learning process through the population and stores 
#' the simulated order of acquisition. As a second step, individuals below a certain threshold are dropped. 
#' The third step is to run OADA using the order of acquisition on a simulated association network 
#' (using the together/apart matrices). When s>0, one can estimate the power of NBDA to 
#' detect social learning for each cut-off point. When s = 0 one can obtain the rate of false positives for each cut-off point.
#' @param x name of the object from the together_apart function 
#' @param sightings group by individual matrix with binary entries 
#' @param cutoff a vector with the values of the desired thresholds 
#' for the inclusion of individuals (number of sightings). 
#' The cut-off point can be as low as 1, while the upper limit has to be 
#' chosen so that at least 2 individuals make the cut-off point.  
#' @param association_index either “SRI” for simple ratio association index or “HWI” for half-weight association index
#' @param iterations how many times the learning process is simulated
#' @param s social learning parameter: if set to 0, one can test for the rate of false positive results (type 1 error). If set to >0, one can test for the power of NBDA.      
#' @param num_ind_learn number of individuals in the population that learn the behaviour. This should be matched with the real number of individuals that learned. 
#' @param cores number of cores the simulations run on in parallel. If unspecified, simulations are set to run on all available cores. 
#' @param keep_learners if FALSE, all individuals are dropped below the specified threshold. If TRUE, all learners are kept regardless of how many times they have been seen. 
#' @param delta_AICc minimum difference in AICc value to select one model (either social or asocial) over the alternative. Defaulted to 2. 
#' @return $raw 	Raw data table with each row representing one model: Columns specify the cut-off point, the number of individuals included, the number of individuals that 
#' acquired the behaviour, the iteration of the respective model, AICc values for both social (aicc) and asocial models (aiccNull), delta AICc, p values, estimates of the 
#' log likelihood of the estimated s and set s, a column specifying if the log likelihood of the estimated s falls within the 95% confidence interval of the set s (yes/no) and 
#' a column specifying if s is an under- or overestimate (under/over) given the log likelihood falls outside the 95% confidence interval. 
#'         $summary	A summary table with models averaged for each cut-off point. Columns specify the cut-off point, the number of individuals included, the number of 
#'         individuals that acquire the behaviour, delta AICc, percentage of models where one model outperforms the other (delta AICc above set threshold), the percentage 
#'         of models where social learning outperforms asocial learning, the set s, averaged estimates and standard deviations for s, proportion of models where the log likelihood 
#'         of s falls within the 95% confidence interval of the set s, a column specifying in what percentage of the cases s is an overestimate or an underestimate respectively 
#'         (given the log likelihood falls outside the 95% confidence interval). 
#' @references \code{citation(Wild and Hoppitt 2019)}

sensitivity_NBDA_ind_error <- function(x, sightings, cutoff, association_index, iterations, s, num_ind_learn, cores=NULL, keep_learners=FALSE, delta_AICc=2){ # define function and parameters
  
  
  
  assoc_SRI <- asnipe::get_network(sightings, data_format="GBI", association_index=association_index) ## calculate association matrix SRI
  
  IDs <- row.names(assoc_SRI) # save names of individuals into a list called 'ID'
  
  
  pb <- tkProgressBar(title = "progress bar", min = 0, # create a progress bar
                      max = length(cutoff), width = 300)
  
  
  #Set number of cores to that in the machine if not specified
  if(is.null(cores)) cores= detectCores(logical = FALSE)
  #Correct number of iterations so it is a multiple of the number of cores
  iterations<-round(iterations/cores)*cores
  
  # split data by ourselves
  chunk.size <- iterations/cores
  
  
  num_sightings <- as.data.frame(colSums(sightings)) ## extract how many times each animal has been seen from the error object
  
  
  
  cutoff_vector <- cutoff # reassign the specified cut off values to a vector with a different name
  
  ## create dataframe to store results
  df7 <- data.frame(matrix(ncol=13, nrow=length(cutoff_vector)*iterations)) ## prepare data frame to store results
  colnames(df7) <- paste(c("cutoff_point", "number_ind", "number_ind_learn", "iteration", "s", "aicc", "aiccNull", "deltaAICc", "p", "loglike_estimated_s", "loglike_set_s", "within_95", "under_over"))
  
  #Save a blank version to be used in parallel processing
  df7Blank<-df7[1:chunk.size,]
  
  # create a dataframe to store results
  results <- data.frame(matrix(nrow=length(cutoff_vector), ncol=12))
  colnames(results) <- paste(c("cutoff", 
                               "num_IDs_total",
                               "num_IDs_learn", 
                               "delta_AICc",
                               "perc_outisde_delta_AICc",
                               "perc_social_learning_supported", 
                               "set_s", 
                               "mean_s", 
                               "sd_s",
                               "within_95%_CI", 
                               "overestimate", 
                               "underestimate"))
  
  for (e in cutoff_vector){ ## run loop for each value of the cutoff vector (e refers to the position of the value in the vector)
    
    
    
    #Set up for parallel processing as detailed in http://www.parallelr.com/r-with-parallel-computing/
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)
    
    res2.p <- foreach(i=1:cores, .combine='rbind') %dopar%
      { 
        
        #re-define the classes and methods used in the loop here      
        source("NBDA code 1.2.15.R")
        
        # local data for results
        res <- df7Blank
        for(m in ((i-1)*chunk.size+1):(i*chunk.size)) {
          
          # create an empty vector where the order of acquisition will be stored in. No length specified
          order_acq <- NULL 
          # create a vector to save the total association with informed individuals
          totalAssocDemon <- rep(length=length(IDs), 0)
          
          # create a vector of length(individuals) filled with 0 to track the status
          statusTracker <- rep(length=length(IDs), 0) 
          
          for (k in 1:num_ind_learn){
            
            
            # fill totalAssocDemon with the association multiplied with the status tracker 
            totalAssocDemon <- assoc_SRI%*%statusTracker
            
            
            Ri <- ((1-statusTracker)*(s*totalAssocDemon+1)) # calculate the learning rate. +1 in the end as no individual-level variables specified
            summedRi <- sum(Ri) # calculate overall learning rate
            learning_probability <- Ri/summedRi # calculate individual learning probability
            
            b <- rmultinom(n=1, size=1, prob=learning_probability) # choose the next individual to learn based on the calculated learning probability
            statusTracker <- as.vector(b)+statusTracker  ## update the status vector with 1 for the individual that learned
            
            order_acq[k] <- which(as.vector(b==1)) # save the individual that learned in a vector order_acq
            
          }
          
          
          if(keep_learners==TRUE){ # if informed individuals should be kept even though they don't make the cutoff
            
            sub <- subset(num_sightings, num_sightings$`colSums(sightings)` >= e) ## drop all animals from list with equal or less than e sightings
            
            list_included <- rownames(sub) ## convert names into a list
            
            OAc1 <- IDs[order_acq] # extract a list of individuals that learned 
            
            list_included_plus <- unique(c(list_included,OAc1))  # create a list of individuals who made the cutoff plus individuals who would not have made it, but acquired the behaviour
            
            OAc2 <- match(OAc1, list_included_plus) # reassing the new positions of these individuals in the reduced social network
            
            IDs_sub <- list_included_plus
            
          } else { # if learners are to be dropped
            
            sub <- subset(num_sightings, num_sightings$`colSums(sightings)` >= e) ## drop all animals from list with equal or less than e sightings
            
            list_included <- rownames(sub) ## convert names into a list
            
            dropped_IDs <- setdiff(IDs, list_included) # extract names of individuals that did not make the cutoff
            dropped_num <- match(dropped_IDs, IDs) # extract the positions of these individuals in the full list of names
            
            
            new_order_acq <- order_acq [! order_acq %in% dropped_num] ## drop individuals from the order of acquistion that did not make the cutoff point
            
            OAc1 <- IDs[new_order_acq] # extract a list of individuals that learned (after dropping the ones that did not make the cutoff)
            OAc2 <- match(OAc1, list_included) # reassing the new positions of these individuals in the reduced social network
            
            IDs_sub <- list_included
          }
          
          assoc_sim <- matrix(rbeta(length(x[,,1]),x[,,1]+1,x[,,2]+1), ncol=length((IDs))) # generate a association matrix based on the number of times dyads have been seen together and seen apart
          # uses a beta distribution B(a+x,b+n-x), where x is the matrix containing the number of times each dyad has been seen together (successes) and n-x a matrix with how many times they have been seen apart (fails)
          # a and b are both 1 (for uniform distribution)
          colnames(assoc_sim) <- IDs
          rownames(assoc_sim) <- IDs
          
          assoc_sim_sub_pre <- assoc_sim[IDs_sub,]
          assoc_sim_sub <- assoc_sim_sub_pre[,IDs_sub]
          
          
          if(length(OAc2!=0)){ # if there are individuals left who learned after removing those below the cutoff
            
            OADA <- oaData(assMatrix=assoc_sim_sub, asoc=NULL, orderAcq=OAc2) ## use the associations with error to run OADA
            OADA_model <- addFit(oadata=OADA, asocialVar=NULL)
            
            param_s <- OADA_model@optimisation$minimum ## corresponds to social learning parameter s
            aicc <- OADA_model@aicc ## extracts aicc of OADA model
            aiccNull <- OADA_model@aiccNull ## extracts aiccNull of OADA model (which is the Nullmodel where s is constrained to 0)
            p <- OADA_model@LRTsocTransPV ## extract p value of the likelihood ratio test
            loglike_s_estimate <- OADA_model@loglik # get log likelihood for the estimated value of s
            loglike_s_set <- addLikelihood(data=OADA, l=s, asocialVar=NULL) # get likelihood for the set s
            
            
            ## save results into local results to be put into dataframe 7
            res[m - (i-1)*chunk.size, "cutoff_point"] <- e
            res[m - (i-1)*chunk.size, "number_ind"] <- length(IDs_sub)
            res[m - (i-1)*chunk.size, "number_ind_learn"] <- length(OAc2)
            res[m - (i-1)*chunk.size,"iteration"] <- m
            res[m - (i-1)*chunk.size,"s"] <- param_s
            res[m - (i-1)*chunk.size,"aicc"] <- aicc
            res[m - (i-1)*chunk.size,"aiccNull"] <- aiccNull
            res[m - (i-1)*chunk.size,"deltaAICc"] <- abs(aicc-aiccNull)
            res[m - (i-1)*chunk.size, "p"] <- p
            res[m - (i-1)*chunk.size,"loglike_estimated_s"] <- loglike_s_estimate
            res[m - (i-1)*chunk.size, "loglike_set_s"] <- loglike_s_set  
            
            # test if the difference between loglikelihood for the set and estimated s is within 1.92 units (95% confidence interval)
            if(abs(loglike_s_set-loglike_s_estimate) <= 1.92){res[m - (i-1)*chunk.size, "within_95"] <- "yes"}  else {
              res[m - (i-1)*chunk.size, "within_95"] <- "no"}
            
            # if it falls outside the 95% confidence interval, calculate the difference between the set s and estimated s to see if the estimated s is an over- or underestimate
            if(res[m - (i-1)*chunk.size,"within_95"]=="yes"){res[m - (i-1)*chunk.size, "under_over"] <- NA} else { 
              if(param_s-s > 0) {res[m - (i-1)*chunk.size,"under_over"] <- "over"} 
              else {res[m - (i-1)*chunk.size,"under_over"] <- "under"}
            }
            
            
            
          }
          
          
          else{
            
            res[m - (i-1)*chunk.size, "cutoff_point"] <- e
            res[m - (i-1)*chunk.size, "number_ind"] <- length(IDs_sub)
            res[m - (i-1)*chunk.size, "number_ind_learn"] <- 0
            res[m - (i-1)*chunk.size,"iteration"] <- m
            res[m - (i-1)*chunk.size,"s"] <- NA
            res[m - (i-1)*chunk.size,"aicc"] <- NA
            res[m - (i-1)*chunk.size,"aiccNull"] <- NA
            res[m - (i-1)*chunk.size,"deltaAICc"] <- NA
            res[m - (i-1)*chunk.size, "p"] <- NA
            res[m - (i-1)*chunk.size,"loglike_estimated_s"] <- NA
            res[m - (i-1)*chunk.size, "loglike_set_s"] <- NA 
            res[m - (i-1)*chunk.size, "within_95"] <- NA 
            res[m - (i-1)*chunk.size, "under_over"] <- NA 
            
            
          }
        }
        res
        
      }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    #Put the results into the df7 object
    df7[(which(cutoff == e)-1)*iterations +1:iterations,]<-res2.p
    
    #temp file that allow us to restart a simulation that crashed partway through
    save(df7, file="sensitivity_NBDA_ind_error_TempOut.dat")
    
    
    # retrieve results for each cut-off point separately
    subset <- subset(df7,df7$cutoff_point==e)
    
    subset2 <- subset(subset, !is.nan(subset[,"deltaAICc"])) # remove rows where AICc difference is NaN (not a number)
    subset2$deltaAICc[subset2[,"deltaAICc"]==Inf] <- 1000 # replace infitinity in the aicc differences with an arbitrary large value
    
    perc_outside_delta_AICc <- (sum(subset2["deltaAICc"]>=delta_AICc, na.rm=TRUE))/iterations*100 # percentage of models where aicc difference is above threshold (excluding NAs)
    subset3 <- subset(subset2,subset2$deltaAICc>=delta_AICc) #create a subset with models only where aicc difference is above the set threshold
    
    perc_support <- (sum(subset3["p"]<0.05, na.rm=TRUE))/iterations*100
    mean_s <- mean(subset3$s, na.rm=TRUE)
    sd_s <- sd(subset3$s, na.rm=TRUE)
    perc_outside <- sum(subset3["within_95"]=="no", na.rm=TRUE )/iterations*100
    
    perc_over <- sum(subset3["under_over"]=="over", na.rm=TRUE)/sum(subset3["within_95"]=="no", na.rm=TRUE)*100
    perc_under <- sum(subset3["under_over"]=="under", na.rm=TRUE)/sum(subset3["within_95"]=="no", na.rm=TRUE)*100
    
    results[which(cutoff==e),1] <- e
    results[which(cutoff==e),2] <- round(mean(subset2$number_ind), digits=2)
    results[which(cutoff==e),3] <- round(mean(subset2$number_ind_learn), digits=2)
    results[which(cutoff==e),4] <- delta_AICc
    results[which(cutoff==e),5] <- perc_outside_delta_AICc
    results[which(cutoff==e),6] <- round(perc_support, digits=2)
    results[which(cutoff==e),7] <- s
    results[which(cutoff==e),8] <- round(mean_s, digits=2)
    results[which(cutoff==e),9] <- round(sd_s, digits=2)
    results[which(cutoff==e),10] <- round(100-perc_outside, digits=2)
    results[which(cutoff==e),11] <- round(perc_over, digits=2)
    results[which(cutoff==e),12] <- round(perc_under, digits=2)
    
    
    setTkProgressBar(pb, which(cutoff_vector==e), label=paste( round(which(cutoff_vector==e)/length(cutoff_vector)*100, 0), "% done"))
  }
  
  
  
  df7 <- df7
  object <- NULL
  object$raw <- df7
  object$summary <- results
  
  return(object)  # return the created object with the two slots $raw and $summary
  close(pb)
  
}


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
