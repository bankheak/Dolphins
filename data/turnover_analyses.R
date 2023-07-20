# Null model for testing population turnover in the humpback dolphins off Richards Bay
# for Shanan Atkins 
# by Mauricio Cantor (m.cantor@ymail.com)
# April 2014, Galapagos Islands

# Summary
# 1. Background
# 2. Package and function loading
# 3. Data input
# 4. Turnover results


# 1. Background -----------------------------------------------------------

# Ultimately, the question here is: Does turnover of dolphins occur during the study period? Broadly, are there units of time where dolphins’ presence is significantly different to that at other, similar units of time? But specifically, this aspect of the analysis is to generate the 95% confidence intervals to give perspective to the (possible) turnover in the Richards Bay humpback dolphin population. This is my interpretation of your question – I trust you will let me know if I haven’t quite got it.

# Null model restrictions
# Dolphins?? sighting frequency for whole8-yr study period (final “grand total” row sum (highlighted yellow)) must stay the same. And when each dolphin was seen within the study period (positions of 1’s across columns and matrices) must be randomised. This is what you called Null Model 1 in your “Disentangling…” paper.


# 2. Package and function loading -----------------------------------------

# loading the package {vegan} which has the function to calculate the Whittaker index

install.packages("vegan", "sfsmisc")

library(vegan)
library(sfsmisc, verbose=F)


# loading the function I've wrote to generate the null distribution and calculate the 95% confidence intervals. Please read the details on the function header and/or open the .html help file attached to this script
# copy the function and paste it in the R console.

#' Population turnover
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

turnover_w=function(data, iter=1000, subseq=FALSE, plot=FALSE){
  
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




# 3. Data input -----------------------------------------------------------

# The raw data: 8 worksheets with 8 sets of binary matrices (of individuals x dates) of my 8-yr study period divided up into
# a. 1 matrix of all 96 months
# b. 2 matrices of 48mo each
# c. 3 matrices of 32mo each
# d. 4 matrices of 24mo each
# e. 6 matrices of 16mo each
# f. 8 matrices of 12mo each
# g. 12 matrices of 8mo each
# h. 16 matrices of 6mo each

# I have created a single matrix for each of 8 sets of binary matrices. In these new matrices, individuals are in the columns and the periods, in the rows. A cell contain the total number of times each individual was sighted within that period. Basically, the rows of the new matrices are the partial column sums (for each block of time) of your previous matrices.

p48m = read.table("p48m.txt",  header= TRUE, row.names=1)
p32m = read.table("p32m.txt",  header= TRUE, row.names=1)
p24m = read.table("p24m.txt",  header= TRUE, row.names=1)
p16m = read.table("p16m.txt",  header= TRUE, row.names=1)
p12m = read.table("p12m.txt",  header= TRUE, row.names=1)
p08m = read.table("p08m.txt",  header= TRUE, row.names=1)
p06m = read.table("p06m.txt",  header= TRUE, row.names=1)
p04m = read.table("p04m.txt",  header= TRUE, row.names=1)
p03m = read.table("p03m.txt",  header= TRUE, row.names=1)

# Transforming into binary matrices
p48m = as.matrix(p48m); p48m[which(p48m>=1)] = 1; p48m[which(p48m<1)] = 0
p32m = as.matrix(p32m); p32m[which(p32m>=1)] = 1; p32m[which(p32m<1)] = 0
p24m = as.matrix(p24m); p24m[which(p24m>=1)] = 1; p24m[which(p24m<1)] = 0
p16m = as.matrix(p16m); p16m[which(p16m>=1)] = 1; p16m[which(p16m<1)] = 0
p12m = as.matrix(p12m); p12m[which(p12m>=1)] = 1; p12m[which(p12m<1)] = 0
p08m = as.matrix(p08m); p08m[which(p08m>=1)] = 1; p08m[which(p08m<1)] = 0
p06m = as.matrix(p06m); p06m[which(p06m>=1)] = 1; p06m[which(p06m<1)] = 0
p04m = as.matrix(p04m); p04m[which(p04m>=1)] = 1; p04m[which(p04m<1)] = 0
p03m = as.matrix(p03m); p03m[which(p03m>=1)] = 1; p03m[which(p03m<1)] = 0



# 4. Turnover results -----------------------------------------------------

t48 = turnover_w(data = p48m, iter = 1000, subseq=F, plot=FALSE)
t32 = turnover_w(data = p32m, iter = 1000, subseq=F, plot=FALSE)
t24 = turnover_w(data = p24m, iter = 1000, subseq=F, plot=FALSE)
t16 = turnover_w(data = p16m, iter = 1000, subseq=F, plot=FALSE)
t12 = turnover_w(data = p12m, iter = 1000, subseq=F, plot=FALSE)
t08 = turnover_w(data = p08m, iter = 1000, subseq=F, plot=FALSE)
t06 = turnover_w(data = p06m, iter = 1000, subseq=F, plot=FALSE)
t04 = turnover_w(data = p04m, iter = 1000, subseq=F, plot=FALSE)
t03 = turnover_w(data = p03m, iter = 1000, subseq=F, plot=FALSE)

all = rbind(t03, t04, t06, t08, t12, t16, t24, t32, t48)
all = cbind(c(03, 04, 06, 08, 12, 16, 24, 32, 48), all)

par(mar=c(4,5,4,1))
# Plot the final results. Whisker represent 95%CI generated by the null model. X-axis represent the number of periods and their respective lengths
errbar(x=c(03, 04, 06, 08, 12, 16, 24, 32, 48), y=all[,2], all[,4], all[,5], ylab="Turnover (Averaged Whittaker Dissimilarity)", pch=1, cap=0.02, xaxt='n', xlab="", las=1, cex=1.0, ylim=c(0.35,0.7), cex.axis=0.8)
axis(1, at=c(03, 04, 06, 08, 12, 16, 24, 32, 48),las=1, cex.axis=0.7)
mtext(side = 1, "Length of periods (months)", line = 2, font = 1)
axis(3, at=c(03, 04, 06, 08, 12, 16, 24, 32, 48),las=1, labels=c(32,24,16,12,8,6,4,3,2), cex.axis=0.7)
mtext(side = 3, "Number of periods", line = 2, font = 1)

# Print final results
all

# I think the population turnover is lower than the null expectation in all cases (Expect for the 48m which had no variance at all and may not make much sense)