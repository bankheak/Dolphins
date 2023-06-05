# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# DYATIC SOCIAL ASSOCIATIONS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

# Load all necessary packages
require(asnipe) # get_group_by_individual--Damien Farine
# Could do permutations
require(assocInd)
require(vegan)
# Run multiple cores for faster computing
require(doParallel)
require(microbenchmark)
require(parallel)
require(foreach)
require(progress)

###########################################################################
# PART 1: Social Association Matrix ---------------------------------------------

# Read in & combine files
firstgen_data <- read.csv("firstgen_data.csv")
secondgen_data <- read.csv("secondgen_data.csv")
orig_data <- rbind(firstgen_data, secondgen_data)

# Make date into a date class
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")
orig_data$Year <- as.numeric(format(orig_data$Date, format = "%Y"))

# Make sure every ID has >10 obs
ID <- unique(orig_data$Code)
obs_vect <- NULL
for (i in 1:length(ID)) {
  obs_vect[i]<- sum(orig_data$Code == ID[i])
}
sub <- data.frame(ID, obs_vect)
sub <- subset(sub, subset=c(sub$obs_vect > 10))
sample_data <- subset(orig_data, orig_data$Code %in% c(sub$ID))
write.csv(sample_data, "sample_data.csv")

# Group each individual by date and sighting
group_data <- cbind(sample_data[,c(2,11,17,21)]) # Seperate date, group and ID
group_data$Group <- cumsum(!duplicated(group_data[1:2])) # Create sequential group # by date
group_data <- cbind(group_data[,3:5]) # Subset ID and group #

# Make a list of only one year per dataframe
years <- unique(group_data$Year)
list_years <- list()
for (i in 1:length(years)) {
  list_years[[i]] <- subset(group_data, subset=c(group_data$Year == years[i]))
}    

# Save nxn list
saveRDS(list_years, file="list_years.RData")

## Test a smaller amount of data for faster results
test <- 100
test_gbi <- get_group_by_individual(list_years[[1]][c(1:100),c(1, 3)], data_format = "individuals")
write.csv(test_gbi, "../data/test_gbi.csv")

# Gambit of the group index
gbi <- list()
for (y in 1:length(years)) {
  gbi[[y]] <- get_group_by_individual(list_years[[y]][,c(1, 3)], data_format = "individuals")
}
saveRDS(gbi, file="gbi.RData")

# Create association matrix
source("../code/functions.R") # SRI & null permutation

n.cores <- detectCores()
system.time({
  registerDoParallel(n.cores)
nxn <- list()
for (i in 1:length(years)) {
  nxn[[i]] <- as.matrix(SRI.func(gbi[[i]]))
}
})

# End parallel processing
stopImplicitCluster()

# Save nxn list
saveRDS(nxn, file="nxn.RData")

###########################################################################
# PART 2: Permutations ---------------------------------------------

# Calculate the CV of the observation association data
# CV = (SD/mean)*100
cv_obs=(sd(nxn) / mean(nxn)) * 100  # Very high CV = unexpectedly high or low association indices in the empirical distribution

#  Create 1000 random group-by-individual binary matrices
reps<- 1000
  registerDoParallel(n.cores)
  nF <- null(gbi, iter=reps)

#' Calculate the association and CV for each of the 1000 permuted matrices to
#' create null distribution
cv_null <- rep(NA,reps)

foreach(i = 1:reps, 
          .combine = c) %dopar% { 
            sri_null = as.matrix(SRI.func(nF[[i]]))
            cv_null[i] <- ( sd(sri_null) / mean(sri_null) ) * 100}
   
stopImplicitCluster()

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
#' from what we would expect by chance. Since the CV(TAI) is lower than the
#' other CV, the associations are lower than expected.
