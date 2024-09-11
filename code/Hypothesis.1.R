# Network Analysis of Anthropogenic Influences in Bottlenose dolphins

# Global Network Analysis Hypothesis #1 #

# Set working directory here
setwd("../data")

# Load all necessary packages
## Predictors
install.packages('igraph', version = '1.6.0') 
library(igraph) # graph_from_adjacency_matrix version = '1.6.0'
if(!require(kinship2)){install.packages('kinship2'); library(kinship2)} # genetic relatedness
if(!require(adehabitatHR)){install.packages('adehabitatHR'); library(adehabitatHR)} # Caluculate MCPs and Kernel density 
## Network
if(!require(ggalt)){install.packages('ggalt'); library(ggalt)} 
if(!require(network)){install.packages('network'); library(network)} # For assigning coordinates to nodes %v%
if(!require(ggmap)){install.packages('ggmap'); library(ggmap)} # register API key version = '3.0.0'
if(!require(graphlayouts)){install.packages('graphlayouts'); library(graphlayouts)} 
if(!require(ggforce)){install.packages('ggforce'); library(ggforce)} # mapping clusters geom_mark_hull
if(!require(ggraph)){install.packages('ggraph'); library(ggraph)} # For network plotting on map
if(!require(tnet)){install.packages('tnet'); library(tnet)} # For weights
if(!require(asnipe)){install.packages('asnipe'); library(asnipe)} # get_group_by_individual
if(!require(assortnet)){install.packages('assortnet'); library(assortnet)} # associative indices
source("../code/functions.R") # nxn
## Mapping
if(!require(statnet)){install.packages('statnet'); library(statnet)}
if(!require(viridis)){install.packages('viridis'); library(viridis)}
if(!require(ggnetwork)){install.packages('ggnetwork'); library(ggnetwork)} # Get cluster coords
if(!require(ggforce)){install.packages('ggforce'); library(ggforce)} # for drawing lines around social clusters
if(!require(ggOceanMaps)){install.packages('ggOceanMaps'); library(ggOceanMaps)} # To map florida
if(!require(intergraph)){install.packages('intergraph'); library(intergraph)} # To use igraph network in ggnet
if(!require(sna)){install.packages('sna'); library(sna)} # For network
if(!require(GGally)){install.packages('GGally'); library(GGally)} # For mapping networks in ggplot version = '2.2.1'
if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)}
if(!require(sf)){install.packages('sf'); library(sf)} # Convert degrees to meters
if(!require(sp)){install.packages('sp'); library(sp)} # Convert degrees to meters
## Bayesian
if(!require(abind)){install.packages('abind'); library(abind)} # array
if(!require(brms)){install.packages('brms'); library(brms)} # For brm model
if(!require(coda)){install.packages('coda'); library(coda)}
if(!require(bayesplot)){install.packages('bayesplot'); library(bayesplot)} # plot parameters in mcmc_area
if(!require(magrittr)){install.packages('magrittr'); library(magrittr)} # For STAN
if(!require(dplyr)){install.packages('dplyr'); library(dplyr)}  # for organizing code
if(!require(rstan)){install.packages('rstan'); library(rstan)} # To make STAN run faster
if(!require(ggrepel)){install.packages('ggrepel'); library(ggrepel)} # for function labs
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(gganimate)){install.packages('gganimate'); library(gganimate)}
if(!require(posterior)){install.packages('posterior'); library(posterior)} # Find the posterior sample names
if(!require(distributional)){install.packages('distributional'); library(distributional)}
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} # Faster computing

# Read in full datasheet and list (after wrangling steps)
orig_data <- read.csv("orig_data.csv") # original data
list_years <- readRDS("list_years.RData") # (1995-2000)/(2001-2006)/(2007-2012)
nxn <- readRDS("nxn.RData") # association matrix of list_years

###########################################################################
# PART 1: Wrangle Data ---------------------------------------------

# Read in & combine files
firstgen_data <- read.csv("firstgen_data.csv") # 1993-2004
secondgen_data <- read.csv("secondgen_data.csv") # 2005-2014
orig_data <- rbind(firstgen_data, secondgen_data) # Combine

# Make date into a date class
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")
orig_data$Year <- as.numeric(format(orig_data$Date, format = "%Y"))

# Fix coding individuals
orig_data$Code <- ifelse(orig_data$Code == "1312", "F222", orig_data$Code)

# Limit to study period
orig_data <- orig_data[orig_data$Year >= 1995 & orig_data$Year <= 2012,]

# Get rid of data without codes
orig_data <- subset(orig_data, subset=c(orig_data$Code != "None"))

# Total Number of individuals in study period
length(unique(orig_data$Code))

# Get rid of any data with no location data
orig_data <- orig_data[!is.na(orig_data$StartLat) & !is.na(orig_data$StartLon),]
orig_data <- subset(orig_data, subset=c(orig_data$StartLat != 999))

# Now split up data 
write.csv(orig_data, "orig_data.csv") # Save data

# Find the codes that appear more than once
repeated_codes <- orig_data$Code[duplicated(orig_data$Code) | duplicated(orig_data$Code, fromLast = TRUE)]

# Get the unique values of these repeated codes
unique_repeated_codes <- unique(repeated_codes)

# Get the length of these unique codes
length(unique_repeated_codes)

## Split data by HAB intensity
orig_data$Group2 <- ifelse(orig_data$Year < 2001, 1, ifelse(orig_data$Year < 2007, 2, 3))
data_by_intense <- split(orig_data, orig_data$Group2) # (1995-2000), (2001-2006), (2007-20012)

# Eliminate IDs with less than 10 locations
fix_list <- function(list_splityears) {
  
  list_years <- list()  # Initialize an empty list to store the updated datasets
  
  for (i in seq_along(list_splityears)) {
    ID <- unique(list_splityears[[i]]$Code)
    obs_vect <- numeric(length(ID))
    
    for (j in seq_along(ID)) {
      obs_vect[j] <- sum(list_splityears[[i]]$Code == ID[j])
    }
    
    sub <- data.frame(ID = ID, obs_vect = obs_vect)
    sub <- subset(sub, subset = obs_vect > 10)
    
    list_years[[i]] <- subset(list_splityears[[i]], Code %in% sub$ID)
    
  }
  
  return(list_years)
}

list_years <- fix_list(data_by_intense)

# Make an overlapping dataset where individuals between period match
overlap_func <- function(list_years) {
## Get unique codes from both lists
codes_list <- list()
for (i in seq_along(list_years)) {
  codes_list[[i]] <- unique(list_years[[i]]$Code)
}

## Find the common codes
common_codes <- Reduce(intersect, codes_list)
## Subset the data frames based on the common codes
list_years_ovrlap <- lapply(list_years, function(df) {
  df[df$Code %in% common_codes, ]
})
return(list_years_ovrlap)
}

list_years <- overlap_func(list_years)

# Save lists
saveRDS(list_years, file = "list_years.RData")

###########################################################################
# PART 2: Social Association Matrix ---------------------------------------------

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
saveRDS(gbi, "gbi.RData")

# Create association matrix
create_nxn <- function(gbi) {
  source("../code/functions.R") # SRI & null permutation
  n.cores <- detectCores()
  system.time({
    registerDoParallel(n.cores)
    nxn <- list()
    for (i in seq_along(gbi)) {
      nxn[[i]] <- as.matrix(SRI.func(gbi[[i]]))
    }                                 
    # End parallel processing
    stopImplicitCluster()
  })
  return(nxn)
}

nxn <- create_nxn(gbi)

# Order data
order_rows <- rownames(nxn[[1]])
order_cols <- colnames(nxn[[1]])

# Apply the order to each matrix in the list
nxn <- lapply(nxn, function(mat) mat[order_rows, order_cols])

# Save nxn lists
saveRDS(nxn, file = "nxn.RData")

# Calculate group size differences
# Read in GBI
gbi <- readRDS("gbi.RData")

# Get the average group size for each ID
group_list <- lapply(gbi, function(group_matrix) {
  
  # Calculate group size for each group
  individual_group_size <- rowSums(group_matrix)
  
  # Create empty vectors to store results
  ids <- character()
  avg_group_sizes <- numeric()
  
  # Iterate through each individual in the group
  for (i in 1:ncol(group_matrix)) {
    
    # Get the individual ID
    individual_id <- colnames(group_matrix)[i]
    
    # Calculate the group size for the individual
    group_size <- ifelse(group_matrix[, individual_id] == 1, 
                         individual_group_size, 0)
    
    # Calculate the average group size for the individual
    avg_group_size <- mean(group_size)
    
    # Append the results to vectors
    ids <- c(ids, individual_id)
    avg_group_sizes <- c(avg_group_sizes, avg_group_size)
  }
  
  # Create a data frame for the current group
  group_data <- data.frame(ID = ids,
                           Average_Group_Size = avg_group_sizes)
  
  return(group_data)
})

# Add HI list
result_df <- readRDS("result_df.RData")
result_df$Group_size <- ifelse(result_df$Period == "1-Before_HAB", 
                               group_list[[1]]$Average_Group_Size[match(result_df$ID, group_list[[1]]$ID)], 
                               ifelse(result_df$Period == "2-During_HAB",
                                      group_list[[2]]$Average_Group_Size[match(result_df$ID, group_list[[2]]$ID)], 
                                      group_list[[3]]$Average_Group_Size[match(result_df$ID, group_list[[3]]$ID)]))

# Change the factor levels and add factor for Period
result_df$Behavior <- ifelse(result_df$HI != "NF", "HI", "None")

# Plot the HI behaviors and group sizes for every year
ggplot(result_df, aes(x = Behavior, y = Group_size, fill = Behavior)) +
  geom_boxplot(outlier.shape = NA) + # Remove outliers
  geom_jitter(aes(color = HI), width = 0.2, alpha = 0.5) + # Set jitter points color and transparency
  facet_wrap(~ Period, labeller = labeller(Period = c("1-Before_HAB" = "Before", 
                                                      "2-During_HAB" = "During", 
                                                      "3-After_HAB" = "After"))) +
  labs(x = "Human-centric Behavior", y = "Individuals' Average Group Size") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank())
# Unequal variance


###########################################################################
# PART 3: CV and Modularity ---------------------------------------------

## Coefficient of Variation ##
# Read in null cv values for one year
cv_null <- readRDS("cv_years.RData")
## Remove NAs, if any
# cv_null = cv_null[!is.na(cv_null)]

# Calculate the CV of the observation association data
# CV = (SD/mean)*100
cv_obs <- lapply(nxn, function (df) {(sd(df) / mean(df)) * 100})  # Very high CV = unexpectedly 
# high or low association indices in the empirical distribution

# Calculate 95% confidence interval, in a two-tailed test
cv_ci = lapply(cv_null, function (df) {quantile(df, probs=c(0.025, 0.975), type=2)})

# Check whether pattern of connections is non-random
par(mfrow=c(3, 1))

# Create a list to store the histograms
hist_cvs <- list()

# Create histograms for each element in cv_null
for (i in seq_along(cv_null)) {
  hist_cvs[[i]] <- hist(cv_null[[i]], 
                        breaks=50,
                        col= "lightblue",
                        main = NULL,
                        xlab="Null CV SRI")
  
  # Add lines for empirical CV, 2.5% CI, and 97.5% CI
  abline(v= cv_obs[[i]], col="red")
  abline(v= cv_ci[[i]], col="blue")
  abline(v= cv_ci[[i]], col="blue")
}

#' This shows whether there are more preferred/avoided 
#' relationships than we would expect at random

## Modularity ##
# Read in data
el <- readRDS("el_years.RData")

## igraph format with weight
n.cores <- detectCores()
system.time({
  registerDoParallel(n.cores)
  dolphin_ig <- list()
  for (j in seq_along(list_years)) {
    dolphin_ig[[j]] <- graph_from_adjacency_matrix(as.matrix(nxn[[j]]),
                                       mode="undirected",
                                       weighted=TRUE, diag=TRUE)
  }  
  ### End parallel processing
  stopImplicitCluster()
})

# Dolphin walk
system.time({
  registerDoParallel(n.cores)
  dolphin_walk <- list()
  for (k in seq_along(dolphin_ig)) {
    dolphin_walk[[k]] <- cluster_walktrap(dolphin_ig[[k]], 
                                          weights = edge_attr(dolphin_ig[[k]], "weight"), 
                                          steps = 4, merges = TRUE, 
                                          modularity = TRUE, 
                                          membership = TRUE)
  } # cluster_edge_betweenness
  
  ### End parallel processing
  stopImplicitCluster()
})

# Run modularity permutations 1000 times for each matrix
run_mod <- function(el, dolphin_walk_list) {
  iter <- 1000
  randmod <- numeric(iter)  # Initialize a numeric vector to store Q-values
  result <- list()
  
  for (k in 1:3) {
    
    for (i in 1:iter) {
    # Save the edgelist into a new object and permutate the link weights
    auxrand <- el[[k]]
    # transform it into igraph format
    igrand <- graph_from_edgelist(auxrand[,1:2]) 
    E(igrand)$weight <- auxrand[,3]
    igrand <- as.undirected(igrand)
    # Now we can permutate the link weights
    E(igrand)$weight <- sample(E(igrand)$weight)
    # calculate the modularity for the permutate copy
    rand_walk <- cluster_walktrap(igrand)
    # and finally save the modularity Q-value into the empty vector
    randmod[i] <- modularity(rand_walk)
  }
  
  # Calculate the 95% confidence interval (two-tailed test)
  ci <- quantile(randmod, probs = c(0.025, 0.975), type = 2)
  
  # Visualization of the random Q distribution
  hist(randmod, xlim = c(0.05, 0.35), main = NA,
       xlab = "Q-value", ylab = "Frequency", col = "lightblue")
  
  # Empirical Q-value
  abline(v = modularity(dolphin_walk_list[[k]]), col = "red")
  
  # 2.5% CI
  abline(v = ci[1], col = "blue")
  
  # 97.5% CI
  abline(v = ci[2], col = "blue")
  
  # Return a data frame with Q-value and confidence intervals
  result[[k]] <- data.frame(Q = modularity(dolphin_walk_list[[k]]), LowCI = ci[1], HighCI = ci[2])
  
  }
  return(result)
}

par(mfrow=c(3, 1))

model <- run_mod(el = el, dolphin_walk_list = dolphin_walk)


###########################################################################
# PART 4: Create ILV and HI Predictors ---------------------------------------------

# SEX and AGE Matrices ------------------------------------------------------

# Read in sex and age data
ILV <- read.csv("Individual_Level_Variables.csv") 
missing_ILV <- read.csv("ILV_missing_age.csv")

# Fix sex so that probable is assigned
ILV$Sex <- ifelse(ILV$Sex == "Probable Female", "Female",
                   ifelse(ILV$Sex == "Probable Male", "Male", ILV$Sex))
 
# Now make a sex and age data frame
ILV_df <- ILV[!duplicated(ILV[, "Alias"]), c("Alias", "Sex", "BirthYear")]
ILV_df$Sex <- ifelse(ILV_df$Sex == "Female", 0, 
                     ifelse(ILV_df$Sex == "Male", 1, NA))
colnames(ILV_df) <- c("Code", "Sex", "Age")
ILV_dem <- ILV_df[ILV_df$Code %in% rownames(nxn[[1]]),]
write.csv(ILV_dem, "ILV_dem.csv")

# Replace missing age
ILV_dem$Age <- ifelse(is.na(ILV_dem$Age), 
                      missing_ILV$Estimated.Birth.Yr..minimum.age.[match(ILV_dem$Code, missing_ILV$Code)], 
                      ILV_dem$Age)
ILV_dem$Age <- ifelse(ILV_dem$Age == "<1988", "1988", ifelse(ILV_dem$Age == "<1997", "1997", ILV_dem$Age))
ILV_dem$Age <- as.numeric(ILV_dem$Age)

# Find the demographics of the population
## Sex
sum(ILV_dem$Sex == "Female")
sum(ILV_dem$Sex == "Male")
## Age
max(ILV_dem$Age)
min(ILV_dem$Age)
hist(ILV_dem$Age)
total_ids <- nrow(ILV_dem)
ids_above_1988 <- sum(ILV_dem$Age >= 1989 & ILV_dem$Age <= 1994)
(ids_above_1988 / total_ids) * 100

# Make sim and diff matrices
ILV_df <- ILV_dem
sim_dif_mat <- function(nxn) {
# Order data
order_rows <- rownames(nxn[[1]])
order_cols <- colnames(nxn[[1]])

# Now reorder the sex and age dataframe
ILV_df <- ILV_df[ILV_df$Code %in% order_rows, ]
ILV_df$Code <- ILV_df$Code[match(order_rows, ILV_df$Code)]

# Make sex sim and age diff matrix
sex_sim <- matrix(0, nrow = nrow(nxn[[1]]), ncol = ncol(nxn[[1]]))
age_diff <- matrix(0, nrow = nrow(nxn[[1]]), ncol = ncol(nxn[[1]]))

for (i in 1:(nrow(sex_sim) - 1)) {
  for (j in (i + 1):ncol(sex_sim)) {
    sex_sim[i, j] <- as.numeric(ILV_df$Sex[i] == ILV_df$Sex[j])
    age_diff[i, j] <- abs(ILV_df$Age[i] - ILV_df$Age[j])
    
    # Since sex_sim and age_diff are symmetric matrices, update the corresponding values
    sex_sim[j, i] <- sex_sim[i, j]
    age_diff[j, i] <- age_diff[i, j]
  }
}
diag(sex_sim) <- 1

return(list(sex_sim, age_diff))
}

ILV_mat <- sim_dif_mat(nxn)

# Normalize dissimilarity to make it range [0,1] and then apply 1-normalized distance
ILV_mat[[2]] <-  1-(ILV_mat[[2]] / max(ILV_mat[[2]]))

# Save ILV matrices
saveRDS(ILV_mat, "ILV_mat.RData")

# HRO Matrix ------------------------------------------------------

# Transform coordinate data into a Spatial Points Dataframe in km
create_coord_data <- function(list, period) {
  
  coords_list <- lapply(list, function(df) {
  
  # Extract IDs and coordinates
  ids <- df$Code
  coordinates <- df[, c("StartLon", "StartLat")]
  
  # Create a SpatialPointsDataFrame with coordinates
  coords_sp <- SpatialPointsDataFrame(coords = coordinates, data = data.frame(id = ids))
  
  # Set CRS to WGS84
  proj4string(coords_sp) <- CRS("+proj=longlat +datum=WGS84")
  
  # Transform to a UTM CRS that uses km as the unit
  coords_sp_utm <- spTransform(coords_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  
  return(coords_sp_utm)})
  
  return(coords_list)
}

dolph.sp <- create_coord_data(list_years)

write.csv(dolph.sp[[1]], "dolph_sp.csv")

# Use the calculated extent in kernelUD
kernel <- lapply(dolph.sp, function(df) kernelUD(df, h = 1000))
saveRDS(kernel, "kernel.RData")
kernel <- readRDS("kernel.RData")

# Calculate Dyadic HRO Matrix: HRO = (Rij/Ri) * (Rij/Rj)
kov <- lapply(kernel, function(df) kerneloverlaphr(df, method = "HR", lev = 95))

# Order data
order_rows <- rownames(nxn[[1]])
order_cols <- colnames(nxn[[1]])

# Apply the order to each matrix in the list
kov <- lapply(kov, function(mat) mat[order_rows, order_cols])

# Save HRO
saveRDS(kov, "kov.RDS")

# GR Matrix ------------------------------------------------------

# Read in sex and age data
ILV_pat <- read.csv("Paternity_data.csv") 

# Order data
order_rows <- rownames(nxn[[1]])
order_cols <- colnames(nxn[[1]])
# Reorder rows in 'ILV' based on 'order_rows'
ILV <- ILV_pat[ILV_pat$Alias %in% order_rows, ]
ILV <- ILV[match(order_rows, ILV$Alias), ]

# Subset paternity data
pedigree_df <- data.frame(Alias = ILV$Alias,
                          Mom = ILV$Mom,
                          Dad = ILV$Dad,
                          Sex = ILV$Sex)

# Fix dad data
pedigree_df$Dad <- ifelse(pedigree_df$Dad == "na", NA, pedigree_df$Dad)
pedigree_df$Dad <- ifelse(pedigree_df$Dad == "FB26 or FB66", "FB26", pedigree_df$Dad)
pedigree_df$Dad <- ifelse(pedigree_df$Dad == "FB76 or FB38", "FB76", pedigree_df$Dad)

# Fix sex so that probable is assigned
pedigree_df$Sex <- ifelse(ILV$Sex == "Probable Female", "Female",
                  ifelse(ILV$Sex == "Probable Male", "Male", ILV$Sex))
# Make sex numeric
pedigree_df$Sex <- ifelse(pedigree_df$Sex == "Female", 2, 
                     ifelse(pedigree_df$Sex == "Male", 1, NA))

# Limit data to non-missing paternity IDs
pedigree_subset <- pedigree_df[!is.na(pedigree_df$Mom) | !is.na(pedigree_df$Dad), ]
# Reset row names to be sequential
row.names(pedigree_subset) <- NULL
saveRDS(pedigree_subset, "pedigree_subset.RData")
pedigree_df <- pedigree_subset

# Make id numeric
## Moms
pedigree_df$ID <- rownames(pedigree_df)
for (i in 1:nrow(pedigree_df)) {
  pedigree_df$Mom <- ifelse(pedigree_df$Mom %in% pedigree_df$Alias[i], 
                            pedigree_df$ID[i], pedigree_df$Mom)
}

## Dads
for (i in 1:nrow(pedigree_df)) {
  pedigree_df$Dad <- ifelse(pedigree_df$Dad %in% pedigree_df$Alias[i], 
                            pedigree_df$ID[i], pedigree_df$Dad)
}

# Only take the ids that aren't found in the 117 list
missing_moms<- subset(pedigree_df, nchar(Mom) > 3)
missing_dads<- subset(pedigree_df, nchar(Dad) > 3)

## Create the sequence of numbers starting from 118
number_mom <- data.frame(Mom = unique(missing_moms$Mom), 
                         ID = c((nrow(pedigree_df) + 1):(nrow(pedigree_df) + length(unique(missing_moms$Mom)))))
## Fill in numbers
for (i in 1:nrow(missing_moms)) {
  missing_moms$Mom <- ifelse(missing_moms$Mom %in% number_mom$Mom[i], 
                             number_mom$ID[i],
                             missing_moms$Mom)
}
## Make ID numeric
missing_moms$Mom <- as.numeric(missing_moms$Mom)

## Do the same thing with dads
number_dad <- data.frame(Dad = unique(missing_dads$Dad), 
                         ID = c((max(missing_moms$Mom) + 1):(max(missing_moms$Mom) + length(unique(missing_dads$Dad)))))
for (i in 1:nrow(missing_dads)) {
  missing_dads$Dad <- ifelse(missing_dads$Dad %in% number_dad$Dad[i], 
                             number_dad$ID[i],
                             missing_dads$Dad)
}
## Make ID numeric
missing_dads$Dad <- as.numeric(missing_dads$Dad)

# Fill in the rest of the NAs with random numbers
## Moms
missing_moms_match <- subset(pedigree_df, nchar(Mom) > 3)
matching_indices <- match(pedigree_df$Mom, missing_moms_match$Mom)
pedigree_df$Mom <- ifelse(!is.na(matching_indices), missing_moms$Mom[matching_indices], pedigree_df$Mom)

## Dads
missing_dads_match<- subset(pedigree_df, nchar(Dad) > 3)
matching_indices <- match(pedigree_df$Dad, missing_dads_match$Dad)
pedigree_df$Dad <- ifelse(!is.na(matching_indices), missing_dads$Dad[matching_indices], pedigree_df$Dad)

# Now create data for function
pedigree_data <- data.frame(id = as.numeric(pedigree_df$ID),
                          mom = as.numeric(pedigree_df$Mom),
                          dad = as.numeric(pedigree_df$Dad),
                          sex = pedigree_df$Sex)

# Assuming your dataframe is named pedigree_data
pedigree_data$dad[is.na(pedigree_data$dad)] <- 0  # Replace NA with 0 or another appropriate code
pedigree_data$mom[is.na(pedigree_data$mom)] <- 0  # Replace NA with 0 or another appropriate code

# Add Fake Fathers
for (i in which(pedigree_data$mom > 0 & pedigree_data$dad == 0)) {
  pedigree_data$dad[i] <- i + max(pedigree_data$dad)
}

# Create fake individuals
fake_ids <- (nrow(pedigree_df) + 1):(max(pedigree_data$dad) + 1)
fake <- data.frame(id = fake_ids,
                   mom = rep(0, length(fake_ids)),
                   dad = rep(0, length(fake_ids)),
                   sex = rep(3, length(fake_ids)))
pedigree_data <- rbind(pedigree_data, fake)

# Change errors
pedigree_data$sex[pedigree_data$id %in% c(139:270)] <- 1
pedigree_data$sex[pedigree_data$id %in% c(118:138)] <- 2

# For limited data
pedigree_data$sex[pedigree_data$id %in% c(94:112, 117:nrow(pedigree_data))] <- 1
pedigree_data$sex[pedigree_data$id %in% c(58:93)] <- 2

# Create GR matrix
ped <- pedigree(id = pedigree_data$id, 
                dadid = pedigree_data$dad, 
                momid = pedigree_data$mom,
                sex = pedigree_data$sex)
plot(ped)

# Calculate kinship matrix
kinship_matrix <- kinship(ped)
kinship_matrix <- kinship_matrix[1:117, 1:117]
saveRDS(kinship_matrix, "kinship_matrix.RData")

# Limited population
kinship_matrix <- kinship_matrix[1:57, 1:57]
saveRDS(kinship_matrix, "kinship_matrix_limit.RData")

# Test correlation to homerange
kinship_matrix <- readRDS("kinship_matrix.RData")
kov <- readRDS("kov.RDS")
mantel(kov[[1]], kinship_matrix)
mantel(kov[[2]], kinship_matrix)
mantel(kov[[3]], kinship_matrix)

# HI Matrices ------------------------------------------------------

# Visualize data: HAB v HI
orig_data$Confirmed_HI <- ifelse(orig_data$ConfHI != "0", 1, 0)

# Filter data where Confirmed_HI == 1
filtered_data <- orig_data[orig_data$Confirmed_HI == 1, ]

# Initialize an empty list to store the counts
yearly_counts <- list()

# Loop through each year
for (year in unique(filtered_data$Year)) {
  # Get the unique individuals for the current year
  unique_ids <- unique(filtered_data$Code[filtered_data$Year == year])
  
  # Store the count of unique IDs in the list
  yearly_counts[[as.character(year)]] <- length(unique_ids)
}

# Convert the list to a data frame for easier viewing
HAB_HI_data <- data.frame(Year = names(yearly_counts), HI_IDs = unlist(yearly_counts))
HAB_HI_data$Year <- as.numeric(HAB_HI_data$Year)
# Add hAB data
HAB_HI_data$HAB <- c(22, 13, rep(0, 2), 5, 0, 12, 8, 18, 5, 38, 19, 2, rep(0, 4), 9)

# Create a barplot
ggplot(aes(x = Year), data = HAB_HI_data) +
  geom_line(aes(y = HI_IDs, color = "HI_IDs", group = 1), size = 1) + 
  geom_point(aes(y = HI_IDs, color = "HI_IDs"), size = 2) + 
  geom_line(aes(y = HAB, color = "HAB", group = 1), size = 1) +  
  geom_point(aes(y = HAB, color = "HAB"), size = 2) + 
  scale_y_continuous(
    name = "Number of human-centric foraging individuals",
    sec.axis = sec_axis(~ ., name = "Number of weeks with >100,000 cells/L")
  ) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010)) +  # Specify the years you want to display
  labs(x = "Year") +
  scale_color_manual(values = c("HI_IDs" = "blue", "HAB" = "orange"), 
                     name = "Variables", 
                     labels = c("Human-centric Behaviors", "Harmful Algal Blooms")) +
  theme(
    panel.background = element_blank()) +
  geom_vline(xintercept = c(2000.5, 2006.5), linetype = "dashed", color = "black", size = 1.5)


# Extract specific columns from each data frame in list_years
aux_data <- function(list_years) {
  aux <- lapply(list_years, function(df) {
    data.frame(Code = df$Code,
               Behaviors = df$Behaviors,
               HumanInteraction = df$HumanInteraction,
               ConfHI = df$ConfHI)})
  
  # Add the 'Foraging' variable to each data frame in the 'aux' list
  aux <- lapply(aux, function(df) {
    df$Foraging <- "Other"
    df$Foraging[grepl(pattern = 'Feed', x = df$Behaviors, ignore.case = FALSE)] <- "Feed"
    df
  })
  return(aux)
}

aux <- aux_data(list_years)

# Categorize ID to Sightings
ID_sight <- function(aux_data) {
  IDbehav <- lapply(aux_data, function(df) {
    ID_by_sightings <- as.data.frame(tapply(df$Code, df$Code, length))
    result <- data.frame(
      Code = rownames(ID_by_sightings),
      Sightings = ID_by_sightings[,1]
    )
    return(result)
  })
  return(IDbehav)
}

IDbehav <- ID_sight(aux)

# Separate HI Behaviors
#' BG = Beg: F, G, H
#' SD = Scavenge and Depredation: A, B, C, D, E
#' FG = Fixed Gear Interaction: P
# Change the code using ifelse statements
subset_HI <- function(aux_data) {
  for (i in seq_along(aux_data)) {
    aux_data[[i]]$DiffHI <- ifelse(aux_data[[i]]$ConfHI %in% c("F", "G", "H"), "BG",
                                   ifelse(aux_data[[i]]$ConfHI %in% c("A", "B", "C", "D", "E"), "SD",
                                          ifelse(aux_data[[i]]$ConfHI %in% c("P"), "FG", 
                                                 ifelse(aux_data[[i]]$Foraging %in% c("Feed")
                                                        & aux_data[[i]]$ConfHI %in% c("0"), "NF", "None"))))
  }
  return(aux_data)  # Return the modified list of data frames
}

aux <- subset_HI(aux)
saveRDS(aux, "aux.RData")
aux <- readRDS("aux.RData")

# Clump all the HI behaviors together
clump_behav <- function(aux_data) {
  for (i in seq_along(aux_data)) {
    aux_data[[i]]$DiffHI <- ifelse(aux_data[[i]]$DiffHI != "None", 1, 0)}
  
  # Categorize DiffHI to IDs
  rawHI <- lapply(aux_data, function(df) {
    
    # Sum up the frequencies of HI by code
    aggregated_df <- aggregate(DiffHI ~ Code, data = df, sum)
    unique_codes_df <- data.frame(Code = unique(df$Code))
    # Merge the unique codes data frame with the aggregated data frame
    merged_df <- merge(unique_codes_df, aggregated_df, by = "Code", all.x = TRUE)
    # Fill missing Freq values (if any) with 0
    merged_df$DiffHI[is.na(merged_df$DiffHI)] <- 0
    
    # Order data
    order_rows <- rownames(nxn[[1]])
    
    # Now reorder the dataframe
    merged_df <- data.frame(Code = order_rows,
                            DiffHI = merged_df$DiffHI[match(order_rows, merged_df$Code)])
    
    return(merged_df)
  })
  return(rawHI)
}

rawHI <- clump_behav(aux)

# Categorize DiffHI to IDs
diff_raw <- function(aux_data) {
  rawHI_diff <- lapply(aux_data, function(df) {
    table_df <- as.data.frame(table(df$Code, df$DiffHI))
    colnames(table_df) <- c("Code", "DiffHI", "Freq")
    
    return(table_df)
  })}

rawHI_diff <- diff_raw(aux)

# Get total number of HI individuals
total_HI_IDs <- unique(unlist(lapply(rawHI, function (df) unique(df$Code[df$DiffHI > 0]))))

# Categorize ID to Sightings
ID_sight <- function(aux_data) {
  IDbehav <- lapply(aux_data, function(df) {
    data <- as.data.frame(table(df$Code))
    colnames(data) <- c("Code", "Sightings")
    # Order data
    order_rows <- rownames(nxn[[1]])
    
    # Now reorder the dataframe
    data <- data %>%
      arrange(match(Code, order_rows))
    
  })
  return(IDbehav)
}

IDbehav <- ID_sight(aux)

# Create a frequency count for each HI behavior
get_IDHI <- function(HI, IDbehav_data, rawHI_diff_data) {
  lapply(seq_along(IDbehav_data), function(i) {
    df <- IDbehav_data[[i]]
    HI_freq <- rawHI_diff_data[[i]]$Freq[rawHI_diff_data[[i]]$DiffHI == HI]
    df$Behav <- HI_freq[match(df$Code, rawHI_diff_data[[i]]$Code)]
    colnames(df) <- c("Code", "Sightings", "Behav")
    return(df)
  })
}

IDbehav_BG <- get_IDHI("BG", IDbehav, rawHI_diff)
IDbehav_FG <- get_IDHI("FG", IDbehav, rawHI_diff)
IDbehav_SD <- get_IDHI("SD", IDbehav, rawHI_diff)
IDbehav_NF <- get_IDHI("NF", IDbehav, rawHI_diff)

saveRDS(IDbehav_BG, "IDbehav_BG.RData")
saveRDS(IDbehav_FG, "IDbehav_BG.RData")
saveRDS(IDbehav_SD, "IDbehav_BG.RData")
saveRDS(IDbehav_NF, "IDbehav_BG.RData")

# Get HI Freq
create_IDbehav_HI <- function(IDbehav_data, rawHI_data){
  IDbehav_HI <- lapply(seq_along(IDbehav_data), function(i) {
    df <- IDbehav_data[[i]]
    df$Behav <- rawHI_data[[i]]$DiffHI
    colnames(df) <- c("Code", "Sightings", "Behav")
    df
  })
  return(IDbehav_HI)
}

IDbehav_HI <- create_IDbehav_HI(IDbehav, rawHI)

# Proportion of Sightings spent in HI
Prop_HI <- function(IDbehav) {
  lapply(seq_along(IDbehav), function(i) {
    df <- IDbehav[[i]]
    df$HIprop <- as.numeric(df$Behav) / as.numeric(df$Sightings)
    df$HIprop[is.na(df$HIprop)] <- 0
    # Keep only 'Code' and 'HIprop' columns
    df <- df[, c('Code', 'HIprop')]
    df
  })
}

prob_HI <- Prop_HI(IDbehav_HI)
prob_BG <- Prop_HI(IDbehav_BG)
prob_SD <- Prop_HI(IDbehav_SD)
prob_FG <- Prop_HI(IDbehav_FG)
prob_NF <- Prop_HI(IDbehav_NF)

saveRDS(prob_BG, "prob_BG.RData")
saveRDS(prob_FG, "prob_FG.RData")
saveRDS(prob_SD, "prob_SD.RData")
saveRDS(prob_NF, "prob_NF.RData")

# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
dis_matr <- function(Prop_HI, nxn) {
  
  # Order data
  order_rows <- rownames(nxn[[1]])
  
  # Apply the order to each matrix in the list
  Prop_HI <- lapply(Prop_HI, function (df) {
    df$Code <- df$Code[match(order_rows, df$Code)]
    return(df)})
  
  # Create matrix
  dissimilarity_HI <- list()
  for (i in seq_along(Prop_HI)) {
    fake_HIprop <- Prop_HI[[i]]$HIprop
    dissimilarity_HI[[i]] <- as.matrix(dist(matrix(fake_HIprop), method = "euclidean"))
  }
  return(dissimilarity_HI)
}

dist_HI <- dis_matr(prob_HI, nxn)

# Transform to similarity
# Method 1: Normalize Euclidian distances to make it range [0,1] and then simply apply 1-normalized distance
sim_HI <- lapply(dist_HI, function (df) {
  # normalize distances to [0,1]
  normdistance <- df / max(df)
  # similarity [0,1] = 1 - distance[0,1]
  similarity1 = 1-normdistance
  return(similarity1)
})

saveRDS(sim_HI, "sim_HI.RData")


###########################################################################
# PART 5: Run MCMC GLMM ---------------------------------------------

# Read in social association matrix and listed data
sim_HI <- readRDS("sim_HI.RData") # HI Sim Matrix
ILV_mat <-readRDS("ILV_mat.RData") # Age and Sex Matrices
kov <- readRDS("kov.RDS")  # Home range overlap
nxn <- readRDS("nxn.RData") # Association Matrix
gr <- readRDS("kinship_matrix.RData")

# Check multicollinearity
mantel(gr, kov[[1]]) # correlated, drop gr

# Prepare random effect for MCMC
num_nodes <- lapply(nxn, function(df) dim(df)[1])
node_names <- lapply(nxn, function(df) colnames(df))

# Separate IDs into i and j
node_ids_i <- lapply(num_nodes, function(df) matrix(rep(1:df, each = df), nrow = df, ncol = df))
node_ids_j <- lapply(node_ids_i, function(df) t(df))

# Format data
upper_tri <- lapply(nxn, function(df) upper.tri(df, diag = TRUE))
edge_nxn <- abind(lapply(nxn, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)

## Split by 3 for int data
HAB_data <- as.data.frame(cbind(c(edge_nxn[,1], edge_nxn[,2], edge_nxn[,3]), 
                                c(rep(1, nrow(edge_nxn)), rep(2, nrow(edge_nxn)), 
                                  rep(3, nrow(edge_nxn)))))
colnames(HAB_data) <- c("SRI", "HAB")
HAB_data$During <- ifelse(HAB_data$HAB == 2, 1, 0)
HAB_data$After <- ifelse(HAB_data$HAB == 3, 1, 0)

HI <- abind(lapply(sim_HI, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
one <- lapply(seq_along(node_ids_i), function(i) factor(as.vector(node_names[[i]][node_ids_i[[i]][upper_tri[[i]]]]), levels = node_names[[i]]))
two <- lapply(seq_along(node_ids_j), function(i) factor(as.vector(node_names[[i]][node_ids_j[[i]][upper_tri[[i]]]]), levels = node_names[[i]]))

# Put data into a dataframe
df_list = data.frame(edge_weight = HAB_data[, 1],
                     HAB_During = HAB_data[, 3],
                     HAB_After = HAB_data[, 4],
                     Period = as.factor(HAB_data[, 2]),
                     HRO = unlist(lapply(kov, function (df) df[upper.tri(df, diag = TRUE)])),
                     sex_similarity = rep(ILV_mat[[1]][upper.tri(ILV_mat[[1]], diag = TRUE)], 3),
                     age_similarity = rep(ILV_mat[[2]][upper.tri(ILV_mat[[2]], diag = TRUE)], 3),
                     #GR = rep(gr[upper.tri(gr, diag = TRUE)], 3),
                     HI_similarity = c(HI[,c(1:3)]),
                     node_id_1 = unlist(one),
                     node_id_2 = unlist(two))

# Make sure that edge_weight is not whole numbers
df_list$edge_weight <- ifelse(df_list$edge_weight == 0, df_list$edge_weight + 0.00001,
                              ifelse(df_list$edge_weight == 1, df_list$edge_weight - 0.00001,
                                     df_list$edge_weight))

# Help STAN run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Multimembership models in brms
fit_sri.0 <- brm(edge_weight ~ (1 | mm(node_id_1, node_id_2)), 
                 family = Beta(), chains = 3, data = df_list)
fit_sri.1 <- brm(edge_weight ~ HI_similarity + Period + 
                   HRO + age_similarity + sex_similarity + 
                   (1 | mm(node_id_1, node_id_2)), 
                 family = Beta(), chains = 3, data = df_list)
fit_sri.2 <- brm(edge_weight ~ HI_similarity * Period +
                   HRO + age_similarity + sex_similarity + 
                   (1 | mm(node_id_1, node_id_2)), 
                 family = Beta(), chains = 3, data = df_list)

# Save data
looic.h1 <- loo(fit_sri.0, fit_sri.1, fit_sri.2, compare = T)
saveRDS(looic.h1, "looic.h1.RData")
looic.h1 <- readRDS("looic.h1.RData")

saveRDS(fit_sri.0, "fit_sri.0.RData")
saveRDS(fit_sri.1, "fit_sri.1.RData")
saveRDS(fit_sri.2, "fit_sri.2.RData")

# Summary Statistics
fit_sri.2 <- readRDS("fit_sri.2.RData")
summary(fit_sri.2)

# Check for model convergence
model <- fit_sri.2
plot(model)
pp_check(model) # check to make sure they line up
# Search how to fix this

# Find the significance
posterior_samples <- as.data.frame(as.matrix( posterior_samples(model) ))
coefficients <- colnames(posterior_samples)
mean(posterior_samples$`b_HI_similarity:Period3` < 0)

# Plot the posterior distribution
# Create mcmc_areas plot
mcmc_plot <- mcmc_intervals(
  as.array(model), 
  pars = c("b_Period3", "b_Period2", 
           "b_HI_similarity:Period3", "b_HI_similarity:Period2",
           "b_HI_similarity",
           "b_HRO", "b_age_similarity", "b_sex_similarity"),
  prob = 0.95, # 95% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) +
  labs(
    title = "Posterior parameter distributions",
    subtitle = "with medians and 95% intervals"
  ) +
  theme_minimal() + # Use a minimal theme
  theme(
    text = element_text(family = "sans"), # Set text family
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(color = "black") # Add axis lines
  )

mcmc_plot + scale_y_discrete(
  labels = c(
    "b_age_similarity" = "Age", 
    "b_sex_similarity" = "Sex", 
    "b_HRO" = "Home-range Overlap", 
    "b_Period2" = "During HAB",
    "b_Period3" = "After HAB",
    "b_HI_similarity" = "Human-centric Similarity",
    "b_HI_similarity:Period2" = "Human-centric Similarity:During HAB",
    "b_HI_similarity:Period3" = "Human-centric Similarity:After HAB"
  )
)


###########################################################################
# PART 6: Display Networks ---------------------------------------------

## Create social network
net <- lapply(nxn, function (df) {
  as.network(df, matrix.type='adjacency',
             directed = F,
             ignore.eval=FALSE,
             names.eval='weight')
})

saveRDS(net, "net.RData")           

# Read in ig object
net <- readRDS("net.RData")
ig <- readRDS("ig.RData")

# Only show IDs of HI dolphins
HI_list <- readRDS("HI_list.RData")
HI_list <- HI_list[-4] # Get rid of natural foragers
HI_IDs <- unique(as.vector(unlist(HI_list))) # Put them all together

#----Modularity---
# igraph format with weight
el_years <- readRDS("el_years.RData")

dolphin_ig <- lapply(nxn, function (mtx) 
  graph_from_adjacency_matrix(as.matrix(mtx),
                              mode="undirected",
                              weighted=TRUE, diag=FALSE))

# Modularity by the WalkTrap algorithm 
dolphin_walk <- lapply(dolphin_ig, function (df)
  cluster_walktrap(df, weights = E(df)$weight,
                   steps = 4, merges = TRUE, 
                   modularity = TRUE, membership = TRUE))

# Create an unweighted network
install.packages('igraph', version = '1.6.0') 
library(igraph)
dolp_ig <- lapply(el_years, function (el) {
  ig <- graph_from_edgelist(el[,1:2])
  # Add the edge weights to this network
  E(ig)$weight <- as.numeric(el[,3])
  # Create undirected network
  ig <- as.undirected(ig)
  return(ig)
  }
)

# Newman's Q modularity
newman <- lapply(dolp_ig, function (df) {
  cluster_leading_eigen(df, steps = -1, weights = E(df)$weight, 
                        start = NULL, options = arpack_defaults(), 
                        callback = NULL, extra = NULL, env = parent.frame())
  })

saveRDS(newman, "newman.RData")
newman <- readRDS("newman.RData")

# Generate a vector of colors based on the number of unique memberships
for (i in seq_along(dolp_ig)) {
  # Generate a vector of colors based on the number of unique memberships
  col <- viridis(min(max(newman[[i]]$membership), length(unique(newman[[i]]$membership))))
  
  # Initialize the color attribute with NA
  V(dolp_ig[[i]])$color <- NA
  
  # Loop through each membership value and assign colors to corresponding vertices
  for (j in 1:max(newman[[i]]$membership)){
    V(dolp_ig[[i]])$color[newman[[i]]$membership == j] <- rep(col[j], sum(newman[[i]]$membership == j))
  }
}

# Read in homerange for individuals
kernel <- readRDS("kernel.RData")

# Create a for loop to store each period's average coordinates
# Extract 50% home range polygons
homerange50 <- lapply(kernel, function (kud) getverticeshr(kud, percent = 50))

# Initialize an empty list to store individual IDs and their centroids
centroid_list <- list()

# Loop through each individual's home range polygons
for(i in seq_along(homerange50)) {
  
  # Initialize an empty data frame to store individual IDs and their centroids
  centroids_df <- data.frame(ID = character(), Latitude = numeric(), Longitude = numeric())
  
  for (id in names(kernel[[1]])) {
    # Convert to sf object for further analysis or export
    homerange50_sf <- st_as_sf(homerange50[[i]], coords = c("X", "Y"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    
    # Get the centroid of the geometry
    centroid <- st_centroid(homerange50_sf$geometry[homerange50_sf$id == id])
    
    # Add the individual ID and centroid coordinates to the data frame
    centroids_df <- rbind(centroids_df, data.frame(ID = id, Latitude = centroid[[1]][2], Longitude = centroid[[1]][1]))
  }
  
  # Put this into a list
  centroid_list[[i]] <- centroids_df
}

# Order data
order_rows <- rownames(nxn[[1]])

centroid_list <- lapply(centroid_list, function(df) { 
  df <- df[df$ID %in% order_rows, , drop = FALSE]  # Subsetting rows based on order_rows
  df <- df[match(order_rows, df$ID), ]  # Reorder rows based on order_rows
  return(df)  # Returning the modified data frame
})

saveRDS(centroid_list, "centroid_list.RData")
centroid_list <- readRDS("centroid_list.RData")

# Define a function to convert UTM coordinates to longitude and latitude
utm_to_lonlat <- function(x, y, zone = 17, northern = TRUE) {
  proj <- sprintf("+proj=utm +zone=%d %s", zone, ifelse(northern, "+north", "+south"))
  xy <- data.frame(x = x, y = y)
  xy <- SpatialPoints(xy, proj4string = CRS(proj))
  xy <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  return(coordinates(xy))
}

# Convert UTM coordinates to longitude and latitude
centroid_list <- lapply(centroid_list, function(df) {
  lonlat <- utm_to_lonlat(df$Longitude, df$Latitude)
  centroid_list <- data.frame(ID = df$ID, 
                              X = lonlat[,1],
                              Y = lonlat[,2])
  return(centroid_list)})

# ---Plot network---
# Set up the plotting area with 1 row and 2 columns for side-by-side plots
labeled_nodes <- list()
plot_list <- list()
register_google(key = "AIzaSyAgFfxIJmkL8LAWE7kHCqSqKBQDvqa9umI")

florida_map <- basemap(limits = c(-87.6349, -79.9743, 24.3963, 31.0006)) + 
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

for (i in 1:length(ig)) {  # Loop through periods
    
    # Load in igraph
    require(igraph)
    
    # Adjust the layout using home range coordinates
    layout_coords <- as.matrix(centroid_list[[i]][, c("X", "Y")])
    adjusted_layout <- layout_coords[order(V(ig[[i]])$name), ]
    
    # Get nodes for each behavior
    labeled_nodes[[i]] <- V(ig[[i]])$name %in% HI_IDs  # Fixed index here

    # Get map of Sarasota, Florida
    sarasota_map <- basemap(limits = c(-82.8, -82.3, 27, 27.6))
    
    # add geographic coordinates
    net_i <- net[[i]]
    net_i %v% "lat" <- layout_coords[,"Y"]
    net_i %v% "lon" <- layout_coords[,"X"]
    x <- net_i %v% "lon"
    y <- net_i %v% "lat"
    
    # Set network and attributes
    node_color <- V(dolp_ig[[i]])$color
    
    # Unrequire igraph
    #detach("package:igraph", unload=TRUE)
    
    # Map the node colors to their corresponding numbers
    color_mapping <- setNames(seq_along(unique(node_color)), unique(node_color))
    node_color_numbers <- color_mapping[node_color]
    grp <- as.vector(node_color_numbers) 
    
    # Graph network
    plot <- ggnetworkmap(
      sarasota_map, # Load in map
      net_i, # Load in network
      size = ifelse(labeled_nodes[[i]], 1.5, 0.5),
      alpha = 0.8, # transparency of nodes
      node.color = ifelse(labeled_nodes[[i]], node_color, "black"), 
      segment.alpha = 0.2, # transparency of edges
      segment.size = get.edge.attribute(net_i, "weight"), # edge thickness
      label.nodes = ifelse(labeled_nodes[[i]], net_i %v% "vertex.names", FALSE),
      label.size = 0.8) +
      geom_encircle(
        aes(x, y, group = as.factor(grp), fill = as.factor(grp)),
        expand = 0.02, 
        alpha = 0.4) +
      scale_fill_brewer(palette = "Set1") +
      theme(
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      ) +
      guides(fill = "none")  # remove legend for fill
    
    plot_list[[i]] <- plot
  
  }

# Plot one at a time
plot_1 <- plot_list[[1]]
pdf("plot_1.pdf", width = 8.5, height = 11)
print(plot_1)
dev.off()
plot_2 <- plot_list[[2]]
pdf("plot_2.pdf", width = 8.5, height = 11)
print(plot_2)
dev.off()
plot_3 <- plot_list[[3]]
pdf("plot_3.pdf", width = 8.5, height = 11)
print(plot_3)
dev.off()

# Create non-mapped graph
labeled_nodes <- list()
plot_list <- list()

for (i in 1:length(ig)) {  # Loop through periods
  
  # Get nodes for each behavior
  labeled_nodes[[i]] <- V(ig[[i]])$name %in% HI_IDs  # Fixed index here
  
  # Make net_i
  net_i <- ig[[i]]
  net_i <- simplify(net_i)
  
  # Set network and attributes
  node_color <- V(dolp_ig[[i]])$color
  
  # Map the node colors to their corresponding numbers
  color_mapping <- setNames(seq_along(unique(node_color)), unique(node_color))
  node_color_numbers <- color_mapping[node_color]
  
  # Update vertex attributes
  V(net_i)$grp <- as.character(node_color_numbers)
  
  create_group_layout <- function(graph, node_groups, layout_func = layout_with_fr, group_spacing = 5) {
  # Get unique groups
  unique_groups <- unique(node_groups)
  
  # Initialize an empty layout
  layout <- matrix(0, nrow = vcount(graph), ncol = 2)
  
  # Apply the layout function to each group individually
  for (i in seq_along(unique_groups)) {
    group_idx <- which(node_groups == unique_groups[i])
    subgraph <- induced_subgraph(graph, group_idx)
    sub_layout <- layout_func(subgraph)
    
    # Normalize sub_layout to avoid large coordinate ranges
    sub_layout <- sub_layout / max(abs(sub_layout))
    
    # Offset sub_layout
    offset_x <- (i - 1) * group_spacing
    offset_y <- (i - 1) * group_spacing
    layout[group_idx, ] <- sub_layout + c(offset_x, offset_y)
  }
  
  return(layout)
}
  
  # Generate a layout based on the node groups
  bb <- create_group_layout(net_i, V(net_i)$grp)
  
  # Create the plot
  plot <- ggraph(net_i, layout = "manual", x = bb[, 1], y = bb[, 2]) +
    geom_edge_link0(aes(color = "grey80", width = E(net_i)$weight), alpha = 0.08) + 
    geom_node_point(aes(fill = grp, size = ifelse(labeled_nodes[[i]], 1.5, 0.5)), shape = 21) +
    geom_text(aes(x = bb[, 1], y = bb[, 2], label = ifelse(labeled_nodes[[i]], V(net_i)$name, "")), color = "black", size = 2, vjust = 1.5) +
    geom_mark_hull(aes(x = bb[, 1], y = bb[, 2], group = grp, fill = grp), 
                   concavity = 4, 
                   expand = unit(2, "mm"), 
                   alpha = 0.25) +
    scale_fill_brewer(palette = "Set1") +
    scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.1), rgb(0, 0, 0, 0.3))) +
    theme_graph() +
    theme(legend.position = "none")
  
  plot_list[[i]] <- plot
}

plot_list[[1]]
plot_list[[2]]
plot_list[[3]]

# What is the cluster size for each period?
combined_cluster_data <- list()

for (i in 1:3) {
  
  ## Get the member data from newman and HI data from labeled_nodes
  member_data <- data.frame(Cluster = newman[[i]]$membership, HI = labeled_nodes[[i]])
  HI_clusters <- member_data[member_data$HI == T,]
  HI_counts <- table(HI_clusters$Cluster)
  membership_counts <- table(newman[[i]]$membership)
  
  ## Ensure both tables have the same keys
  all_keys <- 1:(length(unique(newman[[1]]$membership)))
  
  ## Create a named vector for HI_counts with counts of zero for missing keys
  HI_counts <- as.table(HI_counts)
  HI_counts[as.character(all_keys)] <- HI_counts[as.character(all_keys)]
  HI_counts[is.na(HI_counts)] <- 0
  
  ## Turn data into data frame
  membership_counts_df <- data.frame(Cluster = as.numeric(names(membership_counts)), Total_Cluster_Count = as.numeric(membership_counts))
  HI_counts_df <- data.frame(Cluster = as.numeric(names(HI_counts)), HI_Cluster_Count = as.numeric(HI_counts))
  combined_cluster_data[[i]] <- merge(HI_counts_df, membership_counts_df, all = T)
  combined_cluster_data[[i]]$perc_HI <- (combined_cluster_data[[i]]$HI_Cluster_Count/combined_cluster_data[[i]]$Total_Cluster_Count) *100
  
}

mean(combined_cluster_data[[1]]$Total_Cluster_Count[combined_cluster_data[[1]]$HI_Cluster_Count != 0])
mean(na.omit(combined_cluster_data[[2]]$Total_Cluster_Count[combined_cluster_data[[2]]$HI_Cluster_Count != 0]))
mean(combined_cluster_data[[3]]$Total_Cluster_Count[combined_cluster_data[[3]]$HI_Cluster_Count != 0])

mean(combined_cluster_data[[1]]$HI_Cluster_Count[combined_cluster_data[[1]]$HI_Cluster_Count != 0])
mean(na.omit(combined_cluster_data[[2]]$HI_Cluster_Count[combined_cluster_data[[1]]$HI_Cluster_Count != 0]))
mean(combined_cluster_data[[3]]$HI_Cluster_Count[combined_cluster_data[[1]]$HI_Cluster_Count != 0])

mean(combined_cluster_data[[1]]$perc_HI[combined_cluster_data[[1]]$HI_Cluster_Count != 0])
mean(na.omit(combined_cluster_data[[2]]$perc_HI[combined_cluster_data[[1]]$HI_Cluster_Count != 0]))
mean(combined_cluster_data[[3]]$perc_HI[combined_cluster_data[[1]]$HI_Cluster_Count != 0])
