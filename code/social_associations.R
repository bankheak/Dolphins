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
require(foreach)

###########################################################################
# PART 1: Social Association Matrix ---------------------------------------------

# Read in & combine files
firstgen_data <- read.csv("firstgen_data.csv")
secondgen_data <- read.csv("secondgen_data.csv")
orig_data <- rbind(firstgen_data, secondgen_data)
orig_data <- subset(orig_data, subset=c(orig_data$Code != "None"))

# Make date into a date class
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")
orig_data$Year <- as.numeric(format(orig_data$Date, format = "%Y"))

# Get rid of any data with no location data
sample_data <- orig_data[!is.na(orig_data$StartLat) & !is.na(orig_data$StartLon),]

write.csv(sample_data, "sample_data.csv")
sample_data <- read.csv("sample_data.csv")

# Make a list of three years per dataframe
sample_data$ThreeYearIncrement <- cut(sample_data$Year, breaks = seq(min(sample_data$Year), max(sample_data$Year) + 3, by = 3), labels = FALSE)
list_threeyears <- split(sample_data, sample_data$ThreeYearIncrement)

# Eliminate IDs with less than 5 locations
ID <- list()
for (i in seq_along(list_threeyears)) {
ID[[i]] <- unique(list_threeyears[[i]]$Code)
obs_vect <- NULL
for (j in 1:length(ID[[i]])) {
  obs_vect[j] <- sum(list_threeyears[[i]]$Code == ID[[i]][j])
}
sub <- data.frame(ID = ID[[i]], obs_vect = obs_vect)
sub <- subset(sub, subset=c(sub$obs_vect > 10))
list_threeyears[[i]] <- subset(list_threeyears[[i]], list_threeyears[[i]]$Code %in% c(sub$ID))}

# Save list
saveRDS(list_threeyears, file="list_years.RData")
list_years <- readRDS("list_years.RData")

# Calculate Gambit of the group
gbi <- list()
group_data <- list()
for (i in seq_along(list_years)) {
  
  # Group each individual by date and sighting
  group_data[[i]] <- cbind(list_years[[i]][,c("Date","Sighting","Code","Year")]) 
  group_data[[i]]$Group <- cumsum(!duplicated(group_data[[i]][1:2])) # Create sequential group # by date
  group_data[[i]] <- cbind(group_data[[i]][,3:5]) # Subset ID and group #
  
  # Gambit of the group index
  gbi[[i]] <- get_group_by_individual(group_data[[i]][,c("Code", "Group")], data_format = "individuals")
}

saveRDS(gbi, file="gbi.RData")

# Create association matrix
source("../code/functions.R") # SRI & null permutation

n.cores <- detectCores()
system.time({
  registerDoParallel(n.cores)
nxn <- list()
for (i in seq_along(list_years)) {
  nxn[[i]] <- as.matrix(SRI.func(gbi[[i]]))
}
# End parallel processing
stopImplicitCluster()
})

# Save nxn list
saveRDS(nxn, file="nxn.RData")
nxn <- readRDS("nxn.RData")

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
year = 1
cv_obs=(sd(nxn[[year]]) / mean(nxn[[year]])) * 100  # Very high CV = unexpectedly 
# high or low association indices in the empirical distribution

# Calculate 95% confidence interval, in a two-tailed test
cv_ci = quantile(cv_null[[year]], probs=c(0.025, 0.975), type=2)

# Check whether patterns of connection are non-random
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
#' This shows whether there are more preferred/avoided 
#' relationships than we would expect at random
