# Network Analysis of Anthropogenic Influences in Bottlenose dolphins

# Global Network Analysis Hypothesis #1 #

# Set working directory here
setwd("../data")

# Load all necessary packages
library(asnipe) # get_group_by_individual--Damien Farine
library(assocInd) # Could do permutatioNP
library(vegan)
library(assortnet) # associative indices
library(kinship2) # genetic relatedness
library(ggplot2) # Visualization
library(abind) # array
library(MCMCglmm) # MCMC models
library(coda)
library(bayesplot) # plot parameters
library(sf) # Convert degrees to meters
library(sp) # Creates a SpatialPointsDataFrame by defining the coordinates
library(adehabitatHR) # Caluculate MCPs and Kernel density 
library(rgdal) # Overlap
source("../code/functions.R") # nxn

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

## Split data by HAB intensity
orig_data$Group2 <- ifelse(orig_data$Year < 2001, 1, ifelse(orig_data$Year < 2007, 2, 3))
data_by_intense <- split(orig_data, orig_data$Group2) # (1995-2000), (2001-2006), (2007-20012)

# Eliminate IDs with less than 5 locations
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

###########################################################################
# PART 3: Create ILV and HI Predictors ---------------------------------------------

# SEX and AGE Matrices ------------------------------------------------------

# Read in sex and age data
ILV <- read.csv("Individual_Level_Variables.csv") 

# Fix sex so that probable is assigned
ILV$Sex <- ifelse(ILV$Sex == "Probable Female", "Female",
                   ifelse(ILV$Sex == "Probable Male", "Male", ILV$Sex))
 
# Now make a sex and age data frame
ILV_df <- ILV[!duplicated(ILV[, "Alias"]), c("Alias", "Sex", "BirthYear")]
ILV_df$Sex <- ifelse(ILV_df$Sex == "Female", 0, 
                     ifelse(ILV_df$Sex == "Male", 1, NA))
colnames(ILV_df) <- c("Code", "Sex", "Age")

# Make sim and diff matrices
sim_dif_mat <- function(nxn) {
# Order data
order_rows <- rownames(nxn[[1]])
order_cols <- colnames(nxn[[1]])

# Now reorder the sex and age dataframe
ILV_df <- ILV_df[ILV_df$Code %in% order_rows, ]
ILV_df$Code <- ILV_df$Code[match(order_rows, ILV_df$Code)]

# Estimate unknowns
ILV_df$Sex <- ifelse(is.na(ILV_df$Sex), rbinom(n = nrow(ILV_df), size = 1, prob = 0.5), ILV_df$Sex)
ILV_df$Age <- ifelse(is.na(ILV_df$Age), floor(runif(n = nrow(ILV_df), min = 1970, max = 2002)), as.numeric(ILV_df$Age)) # uniform probability for all ages
ILV_df$Age <- ifelse(is.na(ILV_df$Age), 1997, as.numeric(ILV_df$Age)) # uniform probability for all ages

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

# Use the calculated extent in kernelUD
kernel <- lapply(dolph.sp, function(df) kernelUD(df, h = 1000))

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
ILV <- read.csv("Paternity_data.csv") 

# Order data
order_rows <- rownames(nxn[[1]])
order_cols <- colnames(nxn[[1]])
ILV <- ILV[ILV$Alias %in% order_rows, ]
ILV$Alias <- ILV$Alias[match(order_rows, ILV$Alias)]

# Subset paternity data
pedigree_df <- ILV[!duplicated(ILV[, "Alias"]), c("Alias", "Mom", "Dad", "Sex")]

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

# Make id numeric
pedigree_df$ID <- rownames(pedigree_df)
## Moms
for (i in 1:nrow(pedigree_df)) {
  pedigree_df$Mom <- ifelse(pedigree_df$Mom %in% pedigree_df$Alias[i], 
                            pedigree_df$ID[i], pedigree_df$Mom)
}
### Only take the moms that aren't found in the 117 list
missing_moms<- subset(pedigree_df, nchar(Mom) > 3)
number_mom <- data.frame(Mom = unique(missing_moms$Mom), 
                         ID = c(118:(117 + length(unique(missing_moms$Mom)))))
for (i in 1:nrow(missing_moms)) {
  missing_moms$Mom <- ifelse(missing_moms$Mom %in% number_mom$Mom[i], 
                             number_mom$ID[i],
                             missing_moms$Mom)
}
### Fill in the rest of the NAs with random numbers
missing_moms_match<- subset(pedigree_df, nchar(Mom) > 3)
matching_indices <- match(pedigree_df$Mom, missing_moms_match$Mom)
pedigree_df$Mom <- ifelse(!is.na(matching_indices), missing_moms$Mom[matching_indices], pedigree_df$Mom)

## Dads
for (i in 1:nrow(pedigree_df)) {
  pedigree_df$Dad <- ifelse(pedigree_df$Dad %in% pedigree_df$Alias[i], 
                            pedigree_df$ID[i], pedigree_df$Dad)
}
### Only take the moms that aren't found in the 117 list
missing_dads<- subset(pedigree_df, nchar(Dad) > 3)
number_dad <- data.frame(Dad = unique(missing_dads$Dad), 
                         ID = c(118:(117 + length(unique(missing_dads$Dad)))))
for (i in 1:nrow(missing_dads)) {
  missing_dads$Dad <- ifelse(missing_dads$Dad %in% number_dad$Dad[i], 
                             number_dad$ID[i],
                             missing_dads$Dad)
}
### Fill in the rest of the NAs with random numbers
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

# Get rid of missing individuals in mom and dad
pedigree_data$mom <- ifelse(pedigree_data$mom > 117, 0, pedigree_data$mom)
pedigree_data$dad <- ifelse(pedigree_data$dad > 117, 0, pedigree_data$dad)
pedigree_data$mom <- ifelse(pedigree_data$dad == 0, 0, pedigree_data$mom)
pedigree_data$dad <- ifelse(pedigree_data$mom == 0, 0, pedigree_data$dad)

# Change errors
pedigree_data$sex[pedigree_data$id %in% c(25, 23)] <- 1
pedigree_data$sex[pedigree_data$id %in% c(4, 6, 29)] <- 2

# Create GR matrix
ped <- pedigree(id = pedigree_data$id, 
                dadid = pedigree_data$dad, 
                momid = pedigree_data$mom,
                sex = pedigree_data$sex)
plot(ped)

# Calculate kinship matrix
kinship_matrix <- kinship(ped)

saveRDS(kinship_matrix, "kinship_matrix.RData")

# HI Matrices ------------------------------------------------------

# Visualize data: HAB v HI
HAB_HI_data <- orig_data[, c("Year", "ConfHI")]
HAB_HI_data$ConfHI <- ifelse(HAB_HI_data$ConfHI != "0", 1, 0)
HAB_HI_data <- aggregate(ConfHI ~ Year, data = HAB_HI_data, FUN = function(x) sum(x == 1))
HAB_HI_data$HAB <- c(22, 13, rep(0, 2), 5, 0, 12, 8, 18, 5, 38, 19, 2, rep(0, 4), 9)
# Create a barplot
ggplot(aes(x = Year), data = HAB_HI_data) +
  geom_bar(aes(y = ConfHI, fill = "ConfHI"), stat = "identity", alpha = 0.5, position = position_dodge(width = 0.8)) +
  geom_bar(aes(y = HAB, fill = "HAB"), stat = "identity", alpha = 0.5, position = position_dodge(width = 0.8)) +
  scale_y_continuous(name = "ConfHI", sec.axis = sec_axis(~., name = "HAB")) +
  labs(title = "HAB and HI over Years", x = "Year") +
  scale_fill_manual(values = c("ConfHI" = "blue", "HAB" = "orange"), 
                    name = "Variables", 
                    labels = c("ConfHI", "HAB"))

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
    df <- data.frame(
      Code = unique(df$Code),
      Sightings = tapply(df$Code, df$Code, length)
    )
    df
  })
  return(IDbehav)
}

IDbehav <- ID_sight(aux)

# Separate HI Behaviors
#' BG = Beg: F, G
#' SD = Scavenge and Depredation: A, B, C, D, E
#' FG = Fixed Gear Interaction: P
# Change the code using ifelse statements
subset_HI <- function(aux_data) {
  for (i in seq_along(aux_data)) {
    aux_data[[i]]$DiffHI <- ifelse(aux_data[[i]]$ConfHI %in% c("F", "G"), "BG",
                                   ifelse(aux_data[[i]]$ConfHI %in% c("A", "B", "C", "D", "E"), "SD",
                                          ifelse(aux_data[[i]]$ConfHI %in% c("P"), "FG", "None")))
  }
  return(aux_data)  # Return the modified list of data frames
}

aux <- subset_HI(aux)

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
    return(merged_df)
  })
  return(rawHI)
}

rawHI <- clump_behav(aux)

# Get total number of HI individuals
total_HI_IDs <- unique(unlist(lapply(rawHI, function (df) unique(df$Code[df$ConfHI > 0]))))

# Get HI Freq
create_IDbehav_HI <- function(IDbehav_data, rawHI_data){
  IDbehav_HI <- lapply(seq_along(IDbehav_data), function(i) {
    df <- IDbehav_data[[i]]
    df$HI <- rawHI_data[[i]]$DiffHI
    colnames(df) <- c("Code", "Sightings", "HI")
    df
  })
  return(IDbehav_HI)
}

IDbehav_HI <- create_IDbehav_HI(IDbehav, rawHI)

# Proportion of Sightings spent in HI
Prop_HI <- function(IDbehav) {
  lapply(seq_along(IDbehav), function(i) {
    df <- IDbehav[[i]]
    df$HIprop <- as.numeric(df$HI) / as.numeric(df$Sightings)
    df$HIprop[is.na(df$HIprop)] <- 0
    # Keep only 'Code' and 'HIprop' columns
    df <- df[, c('Code', 'HIprop')]
    df
  })
}

prob_HI <- Prop_HI(IDbehav_HI)

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
# PART 4: Run MCMC GLMM ---------------------------------------------

# Read in social association matrix and listed data
sim_HI <- readRDS("sim_HI.RData") # HI Sim Matrix
ILV_mat <-readRDS("ILV_mat.RData") # Age and Sex Matrices
kov <- readRDS("kov.RDS")  # Home range overlap
nxn <- readRDS("nxn.RData") # Association Matrix
gr <- readRDS("kinship_matrix.RData")

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
HAB_data <- as.data.frame(cbind(c(edge_nxn[,1], edge_nxn[,2], edge_nxn[,3]), c(rep(1, nrow(edge_nxn)), rep(2, nrow(edge_nxn)), rep(3, nrow(edge_nxn)))))
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
                     HRO = unlist(lapply(kov, function (df) df[upper.tri(df, diag = TRUE)])),
                     sex_similarity = rep(ILV_mat[[1]][upper.tri(ILV_mat[[1]], diag = TRUE)], 3),
                     age_similarity = rep(ILV_mat[[2]][upper.tri(ILV_mat[[2]], diag = TRUE)], 3),
                     GR = rep(gr[upper.tri(gr, diag = TRUE)], 3),
                     HI_similarity = c(HI[,c(1:3)]),
                     node_id_1 = unlist(one),
                     node_id_2 = unlist(two))

# Multimembership models in MCMCglmm
fit_mcmc.1 <- MCMCglmm(edge_weight ~ HI_similarity * HAB_During + HI_similarity * HAB_After + 
                         HRO + age_similarity + sex_similarity + GR, 
                       random=~mm(node_id_1 + node_id_2), data = df_list, nitt = 20000) 
summary(fit_mcmc.1)

# Check for model convergence
model <- fit_mcmc.1
plot(model$Sol)
plot(model$VCV)

# Extract Posteriors
posterior <- model$Sol

# Check what proportion is above 0
mean(posterior[, "HI_similarity:HAB_During"] > 0)
mean(posterior[, "HI_similarity:HAB_After"] > 0)

# Plot the posterior distribution
mcmc_intervals(posterior, pars = c("(Intercept)", "HI_similarity", 
                                   "HI_similarity:HAB_During", "HI_similarity:HAB_After",
                                   "HAB_During", "HAB_After", 
                                   "age_similarity", "sex_similarity", "HRO"))
mcmc_areas(
  posterior, 
  pars = c("(Intercept)", 
           "HI_similarity",
           "HI_similarity:HAB_During", 
           "HI_similarity:HAB_After",
           "HAB_During", "HAB_After", 
           "age_similarity", "sex_similarity", "HRO"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

