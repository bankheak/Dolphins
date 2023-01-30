# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# DYATIC SOCIAL ASSOCIATIONS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: Social Association Matrix ---------------------------------------------

# Load all necessary packages
require(asnipe) # get_group_by_individual--Damien Farine
# Could do permutations
require(assocInd)
source("../code/functions.R") # SRI & null permutation

# Read file in
group_data <- read.csv("93-14_data.csv")

# Make date into a date class
group_data$Date <- as.Date(as.character(group_data$Date), format="%d-%b-%y")

# Group each individual by date and sighting
group_data <- cbind(group_data[,c(2,11,17)])
group_data$Group <- cumsum(!duplicated(group_data[1:2]))
group_data <- cbind(group_data[,3:4])
sample_data <- cbind(group_data[,3:4]) # To check that calculations are correct
sample_data<- rbind(sample_data[1:100,])

# Gambit of the group
gbi<- get_group_by_individual(sample_data, data_format = "individuals")

# Create association index
nxn<- SRI.func(gbi)
nxn<-as.matrix(nxn)
write.csv(nxn, "nxn.csv")

# Calculate the CV of the observation association data
# CV = (SD/mean)*100
cv_obs=(sd(nxn) / mean(nxn)) * 100

#  Create 1000 random group-by-individual binary matrices
nF <- null(gbi, iter=1000)

# Calculate the association and CV for each of the 1000 permuted matrices
cv_null <- rep(NA,1000)

for (i in 1:1000) {
  sri_null = as.matrix(SRI.func(nF[[i]]))
  cv_null[i] <- ( sd(sri_null) / mean(sri_null) ) * 100
}

# remove NAs, if any
cv_null = cv_null[!is.na(cv_null)]

# Calculate 95% confidence interval, in a two-tailed test
cv_ci = quantile(cv_null, probs=c(0.025, 0.975), type=2)

# histogram of null CVs
hist(cv_null, 
     breaks=50, 
     col='grey70',
     main = 'Restrictive null model',
     xlab="Null CV SRI")
# empirical CV
abline(v= cv_obs, col="red")
# 2.5% CI
abline(v= cv_ci[1], col="blue")
# 97.5% CI
abline(v= cv_ci[2], col="blue")
#' We can reject the null hypothesis that individuals associate at random
#' and conclude that there is evidence that associations are different 
#' from what we would expect by chance
