# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# DYATIC SOCIAL ASSOCIATIONS
###########################################################################

# Set working directory here
setwd("../data")

# Load all necessary packages
require(asnipe) # get_group_by_individual--Damien Farine
# Could do permutatioNP
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

# Match ID to sex and age data
ILV <- read.csv("Individuals_Residency_Analysis.csv")
## Sex
ID_sex <- setNames(ILV$Sex, ILV$Alias)
orig_data$Sex <- ID_sex[orig_data$Code]
## Age
ID_birth <- setNames(ILV$BirthYear, ILV$Alias)
orig_data$Birth <- ID_birth[orig_data$Code]
orig_data$Age <- as.numeric(orig_data$Year) - as.numeric(orig_data$Birth)

# Get rid of any data with no location data
orig_data <- orig_data[!is.na(orig_data$StartLat) & !is.na(orig_data$StartLon),]
sample_data <- subset(orig_data, subset=c(orig_data$StartLat != 999))

# Get rid of data with no sex or age data
sample_sexage_data <- sample_data[!is.na(sample_data$Sex) & !is.na(sample_data$Age),]

write.csv(sample_data, "sample_data.csv")
sample_data <- read.csv("sample_data.csv")

# Make a list of three years per dataframe
sample_data$ThreeYearIncrement <- cut(sample_data$Year, breaks = seq(min(sample_data$Year), max(sample_data$Year) + 3, by = 3), labels = FALSE)
list_threeyears <- split(sample_data, sample_data$ThreeYearIncrement)

# Make a list of three years per dataframe for sex and age data
sample_sexage_data$ThreeYearIncrement <- cut(sample_sexage_data$Year, breaks = seq(min(sample_sexage_data$Year), max(sample_sexage_data$Year) + 3, by = 3), labels = FALSE)
list_sexage_threeyears <- split(sample_sexage_data, sample_sexage_data$ThreeYearIncrement)

# Eliminate IDs with less than 5 locations
sub_locations <- function(list_years) {
  updated_list_years <- list()  # Initialize an empty list to store the updated datasets
  
  for (i in seq_along(list_years)) {
    ID <- unique(list_years[[i]]$Code)
    obs_vect <- numeric(length(ID))
    
    for (j in seq_along(ID)) {
      obs_vect[j] <- sum(list_years[[i]]$Code == ID[j])
    }
    
    sub <- data.frame(ID = ID, obs_vect = obs_vect)
    sub <- subset(sub, subset = obs_vect > 10)
    
    updated_list_years[[i]] <- subset(list_years[[i]], Code %in% sub$ID)
  }
  return(updated_list_years)
}

list_threeyears <- sub_locations(list_threeyears)
list_sexage_threeyears <- sub_locations(list_sexage_threeyears)

# Save list
saveRDS(list_sexage_threeyears, file="list_sexage_years.RData")
saveRDS(list_threeyears, file="list_years.RData")
list_years <- readRDS("list_years.RData")

# Calculate Gambit of the group
create_gbi <- function(list_years) {
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
          return(gbi)                                      
                                      }

gbi <- create_gbi(list_years)
gbi_sexage <- create_gbi(list_sexage_years)

# Save gbi list
saveRDS(gbi, file="gbi.RData")

# Create association matrix
create_nxn <- function(list_years, gbi) {
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
return(nxn)
}

nxn <- create_nxn(list_years, gbi)
nxn_sexage <- create_nxn(list_sexage_years, gbi_sexage)

# Save nxn lists
saveRDS(nxn, file="nxn.RData")
saveRDS(nxn_sexage, file="nxn_sexage.RData")
nxn <- readRDS("nxn.RData")

###########################################################################
# PART 2: PermutatioNP ---------------------------------------------------------

# Done in the HPC --------------------------------------------------------------

#  Create 1000 random group-by-individual binary matrices
reps<- 1000
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

# Check whether patterNP of connection are non-random
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
#' relatioNPhips than we would expect at random

###########################################################################
# PART 3: SRI Within HI Pairs---------------------------------------------------------

# Read in different behavior's data frames
IDbehav_Beg <- readRDS("IDbehav_Beg.RData")
IDbehav_Pat <- readRDS("IDbehav_Pat.RData")
IDbehav_Dep <- readRDS("IDbehav_Dep.RData")

# Get unique behavior assignments
status <- function(IDbehav, HI, NonHI){
  lapply(seq_along(IDbehav), function(i) {
    IDbehav[[i]]$Stat <- ifelse(IDbehav[[i]]$HI > 0, HI, NonHI)
    df <- IDbehav[[i]][, c('Code', 'Stat')] 
    df
    })
}

## Match each individual with it's behavior
Beg <- status(IDbehav_Beg, "B", "NB")
Pat <- status(IDbehav_Pat, "P", "NP")
Dep <- status(IDbehav_Dep, "D", "ND")

# Replace individuals in the matrix with their assigned behavior
replace_ID_with_HI <- function(sri_matrix, ID_HI_df) {
  # Create vector that matches IDs to their stat
  id_to_stat <- setNames(ID_HI_df$Stat, ID_HI_df$Code)
  
  # Replace each ID with stat in row and column names
  row_names <- id_to_stat[rownames(sri_matrix)]
  col_names <- id_to_stat[colnames(sri_matrix)]
  
  # Create the replaced matrix
  replaced_matrix <- sri_matrix
  
  # Assign row and column names with behavioral states
  dimnames(replaced_matrix) <- list(row_names, col_names)
  return(replaced_matrix)
}

# Make a replaced nxn for each behavior
Beg_nxn <- lapply(seq_along(nxn), function(i) {
  replace_ID_with_HI(nxn[[i]], Beg[[i]])
})
                  
Pat_nxn <- lapply(seq_along(nxn), function(i) {
  replace_ID_with_HI(nxn[[i]], Pat[[i]])
})

Dep_nxn <- lapply(seq_along(nxn), function(i) {
  replace_ID_with_HI(nxn[[i]], Dep[[i]])
})

# Get an average SRI for each category of pairing
## Step 1: Create a matrix for each category of stat

is_NB <- is_B <- list()
for (i in seq_along(Beg_nxn)) {
  is_NB[[i]] <- rownames(Beg_nxn[[i]]) == "NB"
  is_B[[i]] <- rownames(Beg_nxn[[i]]) == "B" 
}

is_NP <- is_P <- list()
for (i in seq_along(Pat_nxn)) {
  is_NP[[i]] <- rownames(Pat_nxn[[i]]) == "NP"
  is_P[[i]] <- rownames(Pat_nxn[[i]]) == "P" 
}

is_ND <- is_D <- list()
for (i in seq_along(Dep_nxn)) {
  is_ND[[i]] <- rownames(Dep_nxn[[i]]) == "ND"
  is_D[[i]] <- rownames(Dep_nxn[[i]]) == "D" 
}

## Step 2: Extract the combinations

### Function to extract combinations
extract_combs <- function(HI_nxn, is_row, is_col) {
  combs <- lapply(seq_along(HI_nxn), function(i) {
    HI_nxn[[i]][is_row[[i]], is_col[[i]]]
  })
  return(combs)
}

#### Apply for each stat comb
NB_NB <- extract_combs(Beg_nxn, is_NB, is_NB)
NB_B <- extract_combs(Beg_nxn, is_NB, is_B)
B_NB <- extract_combs(Beg_nxn, is_B, is_NB)
B_B <- extract_combs(Beg_nxn, is_B, is_B)

NP_NP <- extract_combs(Pat_nxn, is_NP, is_NP)
NP_P <- extract_combs(Pat_nxn, is_NP, is_P)
P_NP <- extract_combs(Pat_nxn, is_P, is_NP)
P_P <- extract_combs(Pat_nxn, is_P, is_P)

ND_ND <- extract_combs(Dep_nxn, is_ND, is_ND)
ND_D <- extract_combs(Dep_nxn, is_ND, is_D)
D_ND <- extract_combs(Dep_nxn, is_D, is_ND)
D_D <- extract_combs(Dep_nxn, is_D, is_D)

## Step 3: Calculate the average of non-diagonal elements in the pairing sub-matrices

### Function to calculate avg
avg_comb <- function(a, b, c, d) {
  avg_a <- lapply(seq_along(a), function(i) {
    mean(a[[i]][lower.tri(a[[i]])])
  })
  avg_b <- lapply(seq_along(b), function(i) {
    mean(b[[i]][lower.tri(b[[i]])])
  })
  avg_c <- lapply(seq_along(c), function(i) {
    mean(c[[i]][lower.tri(c[[i]])])
  })
  avg_d <- lapply(seq_along(d), function(i) {
    mean(d[[i]][lower.tri(d[[i]])])
  })
  avg_df <- data.frame(
    Avg_A = unlist(avg_a),
    Avg_B = unlist(avg_b),
    Avg_C = unlist(avg_c),
    Avg_D = unlist(avg_d)
  )
  return(avg_df)
}


avg_Beg <- avg_comb(NB_NB, NB_B, B_NB, B_B)
colnames(avg_Beg) <- c("NB_NB", "NB_B", "B_NB", "B_B") # Only one beggar in period 4 (2002-2004)
avg_Beg$NB_B <- (avg_Beg$NB_B + avg_Beg$B_NB) / 2
boxplot(avg_Beg[,c(1,2,4)])
plot(avg_Beg[,'B_B'], type="l", col="green", lwd=5, 
     xlab="3-Year Period", ylab="Avg SRI", main = "Beggar-Beggar Pairs")

avg_Pat <- avg_comb(NP_NP, NP_P, P_NP, P_P)
colnames(avg_Pat) <- c("NP_NP", "NP_P", "P_NP", "P_P")
avg_Pat$NP_P <- (avg_Pat$NP_P + avg_Pat$P_NP) / 2
boxplot(avg_Pat[,c(1,2,4)])

avg_Dep <- avg_comb(ND_ND, ND_D, D_ND, D_D)
colnames(avg_Dep) <- c("ND_ND", "ND_D", "D_ND", "D_D") # Only one depredation in period 4 (2002-2004)
avg_Dep$ND_D <- (avg_Dep$ND_D + avg_Dep$D_ND) / 2
boxplot(avg_Dep[,c(1,2,4)])
