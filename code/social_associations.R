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
require(parallel)
require(foreach)

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
sample_data <- read.csv("sample_data.csv")

# Estimate sampling effort and size
## Get estimate of sampling effort
effort <- tapply(sample_data$Date, sample_data$Year, function(x) length(unique(x)))
## Get estimate of population size
unique_ID_year <- tapply(sample_data$Code, sample_data$Year, function(x) length(unique(x)))
## Compare effort to population size
effort <- as.data.frame(effort)
pop <- as.data.frame(unique_ID_year)
pop_effort <- cbind(effort, pop) # Days per year and pop size per year
plot(pop_effort$effort ~ pop_effort$unique_ID_year)
sd(pop_effort$effort)

# Find HI events among individuals
sample_data$ConfHI <- ifelse(sample_data$ConfHI == 0, 0, 1)
ID_HI <- table(sample_data$Code, sample_data$ConfHI, sample_data$Year)
ID_HI <- as.matrix(ID_HI)[,2]

# Group each individual by date and sighting
group_data <- cbind(sample_data[,c("Date","Sighting","Code","Year")]) 
group_data <- subset(group_data, subset=c(group_data$Code != "None"))
group_data$Group <- cumsum(!duplicated(group_data[1:2])) # Create sequential group # by date
group_data <- cbind(group_data[,3:5]) # Subset ID and group #

# Test smaller dataset
test <- 100
group_data <- cbind(group_data[1:test,c(1,3)])
group_data$Year <- rep(c(1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000,
                         2001, 2002),10)

# Make a list of only one year per dataframe
years <- unique(group_data$Year)
list_years <- list()
for (i in 1:length(years)) {
  list_years[[i]] <- subset(group_data, subset=c(group_data$Year == years[i]))
}    

# Save list
saveRDS(list_years, file="list_years.RData")
list_years <- readRDS("list_years.RData")

# Gambit of the group index
gbi <- list()
for (y in 1:length(years)) {
  gbi[[y]] <- get_group_by_individual(list_years[[y]][,c("Code", "Group")], data_format = "individuals")
}
saveRDS(gbi, file="gbi.RData")

##----- Test a smaller amount of data for faster results------
# test_gbi <- gbi[[1]] # or
# test <- 100
# test_gbi <- get_group_by_individual(list_years[[1]][c(1:test),c(1, 3)], data_format = "individuals")
# write.csv(test_gbi, "../data/test_gbi.csv")

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
# PART 2: Permutations ---------------------------------------------------------

# Done in the HPC --------------------------------------------------------------

#  Create 1000 random group-by-individual binary matrices
reps<- 100
  registerDoParallel(n.cores)
  nF <- null(test_gbi, iter=reps)
  saveRDS(nF, "../code/nF.RData")

#' Calculate the association and CV for each of the 1000 permuted matrices to
#' create null distribution
cv_null <- rep(NA,reps)

foreach(i = 1:reps, 
          .combine = c) %dopar% { 
            sri_null = as.matrix(SRI.func(nF[[i]]))
            cv_null[i] <- ( sd(sri_null) / mean(sri_null) ) * 100}
   
stopImplicitCluster()

# Next take results from the HPC ------------------------------------------------

# Read in null cv values for one year
cv_null <- readRDS("../data/cv_years.RData")
## Remove NAs, if any
# cv_null = cv_null[!is.na(cv_null)]

# Calculate the CV of the observation association data
# CV = (SD/mean)*100
cv_obs=(sd(nxn[[year]]) / mean(nxn[[year]])) * 100  # Very high CV = unexpectedly 
# high or low association indices in the empirical distribution

# Calculate 95% confidence interval, in a two-tailed test
cv_ci = quantile(cv_null[[year]], probs=c(0.025, 0.975), type=2)

# histogram of null CVs
hist(cv_null[[year]], 
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
