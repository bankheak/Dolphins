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
# Run multiple cores for faster computing
require(doParallel)
require(microbenchmark)
require(parallel)
require(foreach)
require(progress)

# progress bar ------------------------------------------------------------

iterations <- 100                               # used for the foreach loop  

pb <- progress_bar$new(
  format = "letter = :letter [:bar] :elapsed | eta: :eta",
  total = iterations,    # 100 
  width = 60)

progress_letter <- rep(LETTERS[1:10], 10)  # token reported in progress bar

# allowing progress bar to be used in foreach -----------------------------
progress <- function(n){
  pb$tick(tokens = list(letter = progress_letter[n]))
} 

opts <- list(progress = progress)

###########################################################################
# PART 1: Social Association Matrix ---------------------------------------------

# Read file in
orig_data<- read.csv("secondgen_data.csv")

# Make date into a date class
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")
orig_data$Year <- as.numeric(format(orig_data$Date, format = "%Y"))

# Use only one year
sample_data <- subset(orig_data, subset=c(orig_data$Year == 2005))

# Make sure every ID has >10 obs
ID <- unique(sample_data$Code)
obs_vect <- NULL
for (i in 1:length(ID)) {
  obs_vect[i]<- sum(sample_data$Code == ID[i])
}
sub <- data.frame(ID, obs_vect)
sub <- subset(sub, subset=c(sub$obs_vect > 10))
sample_data <- subset(sample_data, sample_data$Code %in% c(sub$ID))

# Group each individual by date and sighting
group_data <- cbind(sample_data[,c(2,11,17)]) # Seperate date, group and ID
group_data$Group <- cumsum(!duplicated(group_data[1:2])) # Create sequential group # by date
group_data <- cbind(group_data[,3:4]) # Subset ID and group #

# Gambit of the group index
gbi<- get_group_by_individual(sample_data, data_format = "individuals")
write.csv(gbi, "gbi.csv")

# Create association matrix
source("../code/functions.R") # SRI & null permutation

n.cores <- detectCores()
system.time({
  registerDoParallel(n.cores)
  nxn<- SRI.func(gbi)
})
nxn<-as.matrix(nxn)
# end parallel processing
stopImplicitCluster()


###########################################################################
# PART 2: Permutations ---------------------------------------------

# Calculate the CV of the observation association data
# CV = (SD/mean)*100
cv_obs=(sd(nxn) / mean(nxn)) * 100  # Very high CV = unexpectedly high or low association indices in the empirical distribution

#  Create 1000 random group-by-individual binary matrices
reps<- 1000
nF <- null(gbi, iter=reps)

#' Calculate the association and CV for each of the 1000 permuted matrices to
#' create null distribution
cv_null <- rep(NA,reps)

##' try to put in parallel
##' trial<- lapply(nF, SRI.func)
##' sri_null = as.matrix(trial)
##' look-up foreach progess bar ##
##' Look at cyberduck sstp to run faster computing 2 ##


registerDoParallel(n.cores)

  null <- foreach(i = 1:reps, 
          .combine = c) %dopar% { 
            # replace c with rbind to create a dataframe
            sri_null = as.matrix(SRI.func(nF[[i]]))
            ( sd(sri_null) / mean(sri_null) ) * 100}
   
stopImplicitCluster()

for (i in 1:reps) {
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
#' from what we would expect by chance. Since the CV(TAI) is lower than the
#' other CV, the associations are lower than expected.
