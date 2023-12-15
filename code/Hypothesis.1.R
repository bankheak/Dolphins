# Network Analysis of Anthropogenic Influences in Bottlenose dolphins

# Set working directory here
setwd("../data")

# Load all necessary packages
require(asnipe) # get_group_by_individual--Damien Farine
require(assocInd) # Could do permutatioNP
require(vegan)
require(doParallel) # Run multiple cores for faster computing
require(foreach)
library(reshape2) # For graphing
library(gridExtra) # To combine plots

# Read in full datasheet and list (after wrangling steps)
orig_data <- read.csv("orig_data.csv")
list_years <- readRDS("list_years_ovrlap.RData")
nxn <- readRDS("nxn_ovrlap.RData")

###########################################################################
# PART 1: Wrangle Data ---------------------------------------------

# Read in & combine files
firstgen_data <- read.csv("firstgen_data.csv") # 1993-2004
secondgen_data <- read.csv("secondgen_data.csv") # 2005-2014
orig_data <- rbind(firstgen_data, secondgen_data) # Combine

# Get rid of data without codes
orig_data <- subset(orig_data, subset=c(orig_data$Code != "None"))

# Get rid of any data with no location data
orig_data <- orig_data[!is.na(orig_data$StartLat) & !is.na(orig_data$StartLon),]
orig_data <- subset(orig_data, subset=c(orig_data$StartLat != 999))

# Make date into a date class
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")
orig_data$Year <- as.numeric(format(orig_data$Date, format = "%Y"))

# Fix coding individuals
orig_data$Code <- ifelse(orig_data$Code == "1312", "F222", orig_data$Code)

# Now split up data 7 years before and after HAB
orig_data <- orig_data[orig_data$Year >= 1998 & orig_data$Year <= 2011,]
write.csv(orig_data, "orig_data.csv") # Save data

# Make a list of split years per dataframe for before and after HAB
## Sort the data by Year
orig_data <- orig_data[order(orig_data$Year), ]
## Create a column indicating the group (1 or 2) based on the midpoint
orig_data$Group <- ifelse(orig_data$Year < median(orig_data$Year), 1, 2)
## Split the data into two groups based on the 'Group' column
list_splityears <- split(orig_data, orig_data$Group)
## Split data by year
data_by_year <- split(orig_data, orig_data$Year)

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
  # Fix sex so that probable is assigned
  list_years <- lapply(list_years, function(df) {
    df$Sex <- ifelse(df$Sex == "Probable Female", "Female",
                     ifelse(df$Sex == "Probable Male", "Male", df$Sex))
    return(df)
  })
  
  return(list_years)
}

list_years <- fix_list(list_splityears)
list_years_one <- fix_list(data_by_year)

# Make an overlapping dataset where individuals before and after HAB match
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

## Subset the data_by_year based on the common codes
list_years_one_ovrlap <- lapply(list_years_one, function(df) {
  df[df$Code %in% unique(list_years_ovrlap[[1]]$Code), ]
})

# Save lists
saveRDS(list_years_ovrlap, file = "list_years_ovrlap.RData")
saveRDS(list_years_one_ovrlap, file = "list_years_one_ovrlap.RData")


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

gbi <- create_gbi(list_years_one_ovrlap)
gbi_ovrlap <- create_gbi(list_years_ovrlap)

# Save gbi lists
saveRDS(gbi, file = "gbi.RData")
saveRDS(gbi_ovrlap, file = "gbi_ovrlap.RData")

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
nxn_ovrlap <- create_nxn(gbi_ovrlap)

# Save nxn lists
saveRDS(nxn, file = "nxn.RData")
saveRDS(nxn_ovrlap, file = "nxn_ovrlap.RData")

# Calculate nxn SE over years
calc_cell_se <- function(matrices_list, se_data) {
  
  # Get number of years per matrix
  num_year <- length(se_data)
  
  # Get the dimensions of the matrices
  mat_dims <- dim(matrices_list)
  
  # Codes
  codes <- colnames(matrices_list)
  
  # Initialize a matrix to store standard errors
  se_matrix <- matrix(0, nrow = mat_dims[1], ncol = mat_dims[2])
  
  # Calculate standard error for each cell
  for (i in 1:mat_dims[1]) {
    for (j in 1:mat_dims[2]) {
      cell_data <- lapply(se_data, function(matrix) matrix[rownames(matrix)[i,] %in% codes, colnames(matrix[,j]) %in% codes])
      se_matrix[i, j] <- sd(cell_data) / sqrt(num_year)
    }
  }
  
  return(se_matrix)
}

# Calculate standard error for each cell in the matrices
se_matrix1 <- calc_cell_se(matrices_list = nxn_ovrlap[[1]], se_data = list(nxn[c(1:7)]))
se_matrix2 <- calc_cell_se(matrices_list = nxn_ovrlap[[2]], se_data = list(nxn[c(8:14)]))
SE_array <- list(SE_array1, SE_array2)

saveRDS(SE_array, file = "SE_array.RData")

###########################################################################
# PART 3: Create ILV and HI Predictors ---------------------------------------------

# SEX and AGE Matrices ------------------------------------------------------

# Read in sex and age data
ILV <- read.csv("Individual_Level_Variables.csv") 

# Now make a sex and age data frame
ILV_df <- list_years[[2]][!duplicated(list_years[[2]][, "Code"]), c("Code", "Sex", "Age")]
ILV_df$Sex <- ifelse(ILV_df$Sex == "Female", 0, 
                     ifelse(ILV_df$Sex == "Male", 1, NA))

# Order data
order_rows <- rownames(nxn[[1]])
order_cols <- colnames(nxn[[1]])

# Apply the order to each matrix in the list
nxn <- lapply(nxn, function(mat) mat[order_rows, order_cols])
SE_list <- lapply(SE_list, function(mat) mat[order_rows, order_cols])

# Now reorder the sex and age dataframe
ILV_df$Code <- ILV_df$Code[match(order_rows, ILV_df$Code)]

# Estimate unknowns
ILV_df$Sex <- ifelse(is.na(ILV_df$Sex), rbinom(n = nrow(ILV_df), size = 1, prob = 0.5), ILV_df$Sex)
ILV_df$Age <- ifelse(is.na(ILV_df$Age), floor(runif(n = nrow(ILV_df), min = 1, max = 40)), ILV_df$Age) # uniform probability for all ages

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

# HRO Matrix ------------------------------------------------------

# Read in file
list_years <- readRDS("list_years_ovrlap.RData")

# Aggregate list into one homerange overlap matrix
list_years_df <- merge(list_years[[1]], list_years[[2]], all = T)

# Transform coordinate data into a Spatial Points Dataframe in km
create_coord_data <- function(df) {
  
  # Extract IDs and coordinates
  ids <- df$Code
  coordinates <- df[, c("StartLon", "StartLat")]
  
  # Create a SpatialPointsDataFrame with coordinates
  coords_sp <- SpatialPointsDataFrame(coords = coordinates, data = data.frame(id = ids))
  
  # Set CRS to WGS84
  proj4string(coords_sp) <- CRS("+proj=longlat +datum=WGS84")
  
  # Transform to a UTM CRS that uses km as the unit
  coords_sp_utm <- spTransform(coords_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  
  return(coords_sp_utm)
}

dolph.sp <- create_coord_data(list_years_df)

# Visualize data extent
vis.sf <- function(dolph.sp) {
  dolp.sf <- st_as_sf(dolph.sp)
  return(dolp.sf)
}

dolph.sf <- vis.sf(dolph.sp)

vis_coord <- lapply(seq_along(dolph.sf), function (i) {ggplot(dolph.sf[[i]]) +
    geom_sf(aes(color = "Dolphins"), size = 2, alpha = 0.5) +
    theme_bw() +
    labs(title = "Distribution") +
    scale_color_manual(values = c("Dolphins" = "blue"))})
wrap_plots(vis_coord, nrow = 1)

# Look into what bandwidth 
value <- 1
period <- 1
n <- length(dolph.sp[[period]]@coords[, value]) 
bw_scott <- (4 / (3 * n))^(1/5) * sd(dolph.sp[[period]]@coords[, value])  # Scott's rule
## Use the calculated bandwidth in density estimation
density_estimate <- density(dolph.sp[[period]]@coords[, value], bw = bw_scott)
## Plot the density estimate
plot(density_estimate)

# Calculate the maximum distance between points
max(sqrt((dolph.sp[[period]]@coords[,1] - min(dolph.sp[[period]]@coords[,1]))^2 + 
           (dolph.sp[[period]]@coords[,2] - min(dolph.sp[[period]]@coords[,2]))^2))

# Use the calculated extent in kernelUD
kernel <- kernelUD(dolph.sp, h = 1000)

# Calculate Dyadic HRO Matrix: HRO = (Rij/Ri) * (Rij/Rj)
kov <- kerneloverlaphr(kernel, method = "HR", lev = 95)

# Save HRO
saveRDS(kov, "kov.RDS")
# GR Matrix ------------------------------------------------------
# HI Matrices ------------------------------------------------------

# Separate HI Behaviors
#' BG = Beg: F, G
#' SD = Scavenge and Depredation: B, C, D, E, A
#' FG = Fixed Gear Interaction: P
# Change the code using ifelse statements
subset_HI <- function(aux_data) {
  for (i in seq_along(aux_data)) {
    aux_data[[i]]$DiffHI <- ifelse(aux_data[[i]]$ConfHI %in% c("F", "G"), "BG",
                                   ifelse(aux_data[[i]]$ConfHI %in% c("B", "C", "D", "E", "A"), "SD",
                                          ifelse(aux_data[[i]]$ConfHI %in% c("P"), "FG", "None")))
  }
  return(aux_data)  # Return the modified list of data frames
}

aux <- subset_HI(aux)

# Categorize DiffHI to IDs
diff_raw <- function(aux_data) {
  rawHI_diff <- lapply(aux_data, function(df) {
    table_df <- as.data.frame(table(df$Code, df$DiffHI))
    colnames(table_df) <- c("Code", "DiffHI", "Freq")
    return(table_df)
  })}

rawHI_diff <- diff_raw(aux)

# Create a frequency count for each HI behavior
get_IDHI <- function(HI, IDbehav_data, rawHI_diff_data) {
  lapply(seq_along(IDbehav_data), function(i) {
    df <- IDbehav_data[[i]]
    HI_freq <- rawHI_diff_data[[i]]$Freq[rawHI_diff_data[[i]]$DiffHI == HI]
    df$HI <- HI_freq[match(df$Code, rawHI_diff_data[[i]]$Code)]
    colnames(df) <- c("Code", "Foraging", "HI")
    df
  })
}

# Including zeros
IDbehav_BG <- get_IDHI("BG", IDbehav, rawHI_diff)
IDbehav_FG <- get_IDHI("FG", IDbehav, rawHI_diff)
IDbehav_SD <- get_IDHI("SD", IDbehav, rawHI_diff)

# Clump all the HI behaviors together
clump_behav <- function(aux_data) {
  for (i in seq_along(aux_data)) {
    aux_data[[i]]$ConfHI <- ifelse(aux_data[[i]]$ConfHI != "0", 1, 0)}
  
  # Categorize ConfHI to IDs
  rawHI <- lapply(aux_data, function(df) {
    # Sum up the frequencies of HI by code
    aggregated_df <- aggregate(ConfHI ~ Code, data = df, sum)
    unique_codes_df <- data.frame(Code = unique(df$Code))
    # Merge the unique codes data frame with the aggregated data frame
    merged_df <- merge(unique_codes_df, aggregated_df, by = "Code", all.x = TRUE)
    # Fill missing Freq values (if any) with 0
    merged_df$ConfHI[is.na(merged_df$ConfHI)] <- 0
    return(merged_df)
  })
  return(rawHI)
}

rawHI <- clump_behav(aux)

# Get HI Freq
create_IDbehav_HI <- function(IDbehav_data, rawHI_data){
  IDbehav_HI <- lapply(seq_along(IDbehav_data), function(i) {
    df <- IDbehav_data[[i]]
    df$HI <- rawHI_data[[i]]$ConfHI
    colnames(df) <- c("Code", "Foraging", "HI")
    df
  })
  return(IDbehav_HI)
}

IDbehav_HI <- create_IDbehav_HI(IDbehav, rawHI)

# Proportion of time Foraging spent in HI
Prop_HI <- function(IDbehav) {
  lapply(seq_along(IDbehav), function(i) {
    df <- IDbehav[[i]]
    df$HIprop <- as.numeric(df$HI) / as.numeric(df$Foraging)
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

# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
dis_matr <- function(Prop_HI) {
  dissimilarity_HI <- list()
  for (i in seq_along(Prop_HI)) {
    fake_HIprop <- Prop_HI[[i]]$HIprop
    dissimilarity_HI[[i]] <- as.matrix(dist(matrix(fake_HIprop), method = "euclidean"))
  }
  return(dissimilarity_HI)
}

dist_HI <- dis_matr(prob_HI)
dist_BG <- dis_matr(prob_BG)
dist_SD <- dis_matr(prob_SD)
dist_FG <- dis_matr(prob_FG)

saveRDS(dist_HI, "dist_HI.RData")
saveRDS(dist_BG, "dist_BG.RData")
saveRDS(dist_SD, "dist_SD.RData")
saveRDS(dist_FG, "dist_FG.RData")




###########################################################################
# PART 4: Run MCMC GLMM ---------------------------------------------

# Prepare random effect for MCMC
num_nodes <- lapply(nxn, function(df) dim(df)[1])
node_names <- lapply(nxn, function(df) colnames(df))

# Separate IDs into i and j
node_ids_i <- lapply(num_nodes, function(df) matrix(rep(1:df, each = df), nrow = df, ncol = df))
node_ids_j <- lapply(node_ids_i, function(df) t(df))

# Format data
upper_tri <- lapply(nxn, function(df) upper.tri(df, diag = TRUE))
edge_nxn <- abind(lapply(nxn, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
HAB_data <- cbind(c(edge_nxn[,1], edge_nxn[,2]), c(rep(0, nrow(edge_nxn)), rep(1, nrow(edge_nxn))))
HI <- abind(lapply(dist_HI, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
BG <- abind(lapply(dist_BG, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
FG <- abind(lapply(dist_FG, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
SD <- abind(lapply(dist_SD, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
#SE <- abind(lapply(SE_list, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
one <- lapply(seq_along(node_ids_i), function(i) factor(as.vector(node_names[[i]][node_ids_i[[i]][upper_tri[[i]]]]), levels = node_names[[i]]))
two <- lapply(seq_along(node_ids_j), function(i) factor(as.vector(node_names[[i]][node_ids_j[[i]][upper_tri[[i]]]]), levels = node_names[[i]]))

# Put data into a dataframe
df_list = data.frame(edge_weight = HAB_data[, 1],
                     HAB = HAB_data[, 2],
                     HRO = rep(kov[upper.tri(kov, diag = TRUE)], 2),
                     sex_similarity = rep(sex_sim[upper.tri(sex_sim, diag = TRUE)], 2),
                     age_difference = rep(scale(c(age_diff[upper.tri(age_diff, diag = TRUE)])), 2),
                     #GR = gr_list,
                     HI_differences = c(HI[,1], HI[,2]),
                     BG_differences = c(BG[,1], BG[,2]),
                     FG_differences = c(FG[,1], FG[,2]),
                     SD_differences = c(SD[,1], SD[,2]),
                     #Obs.Err = c(SE[,1], SE[,2]),
                     node_id_1 = c(one[[1]], one[[2]]),
                     node_id_2 = c(two[[1]], two[[2]]))

# Multimembership models in MCMCglmm
## HI Behavior Combined
fit_mcmc.1 <- MCMCglmm(edge_weight ~ HI_differences * HAB + HRO + age_difference + sex_similarity, 
                     random=~mm(node_id_1 + node_id_2), data = df_list, nitt = 20000)
summary(fit_mcmc)

## HI Behavior Separated
fit_mcmc.2 <- MCMCglmm(edge_weight ~ HI_differences * HAB + HRO + age_difference + sex_similarity, 
                     random=~mm(node_id_1 + node_id_2), data = df_list, nitt = 20000)
summary(fit_mcmc)

# Check for model convergence
plot(fit_mcmc$Sol)
plot(fit_mcmc$VCV)

# Extract Posteriors
posterior <- fit_mcmc$Sol

# Plot the posterior distribution
mcmc_intervals(posterior, pars = c("(Intercept)", "HI_differences", "HI_differences:HAB", "HAB", 
                                   "age_difference", "sex_similarity"))
mcmc_areas(
  posterior, 
  pars = c("(Intercept)", "HI_differences:HAB", "HAB", 
           "age_difference", "sex_similarity"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

