# Network Analysis of Anthropogenic Influences in Bottlenose dolphins

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
list_years <- readRDS("list_years.RData") # (1998-2004)/(2005-2014)
list_years_int <- readRDS("list_years_int.RData") # (1995-2000)/(2001-2006)/(2007-20012)
nxn <- readRDS("nxn.RData") # association matrix of list_years
nxn_int <- readRDS("nxn_int.RData") # association matrix of list_years_int

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

# Now split up data 
orig_data <- orig_data[orig_data$Year >= 1995 & orig_data$Year <= 2012,]
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

list_years <- fix_list(list_splityears)
list_years_one <- fix_list(data_by_year)
list_years_int <- fix_list(data_by_intense)

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
list_years_int <- overlap_func(list_years_int)

## Subset the data_by_year based on the common codes
SE_split <- function(list_years) {
list_years_one <- lapply(list_years_one, function(df) {
  df[df$Code %in% unique(list_years[[1]]$Code), ]
})}

list_years_one <- SE_split(list_years)
list_years_one_int <- SE_split(list_years_int)

# Save lists
saveRDS(list_years, file = "list_years.RData")
saveRDS(list_years_int, file = "list_years_int.RData")
saveRDS(list_years_one, file = "list_years_one.RData")
saveRDS(list_years_one, file = "list_years_one_int.RData")


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
gbi_one <- create_gbi(list_years_one)
gbi_int <- create_gbi(list_years_int)

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
nxn_one <- create_nxn(gbi_one)
nxn_int <- create_nxn(gbi_int)

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
## Two
se_matrix1 <- calc_cell_se(matrices_list = nxn[[1]], se_data = list(nxn_one[c(1:7)]))
se_matrix2 <- calc_cell_se(matrices_list = nxn[[2]], se_data = list(nxn_one[c(8:14)]))
SE_list <- list(se_matrix1, se_matrix2)
## Three
se_matrix1 <- calc_cell_se(matrices_list = nxn[[1]], se_data = list(nxn_one[c(1:7)]))
se_matrix2 <- calc_cell_se(matrices_list = nxn[[2]], se_data = list(nxn_one[c(8:14)]))
SE_list <- list(se_matrix1, se_matrix2)

# Order data
order_rows <- rownames(nxn[[1]])
order_cols <- colnames(nxn[[1]])

# Apply the order to each matrix in the list
nxn <- lapply(nxn, function(mat) mat[order_rows, order_cols])
SE_list <- lapply(SE_list, function(mat) mat[order_rows, order_cols])

# Save nxn lists
saveRDS(nxn, file = "nxn.RData")
saveRDS(nxn_one, file = "nxn_one.RData")
saveRDS(nxn_int, file = "nxn_int.RData")

saveRDS(SE_list, file = "SE_list.RData")
saveRDS(SE_list_int, file = "SE_list_int.RData")

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
ILV_mat_int <- sim_dif_mat(nxn_int)

# Save ILV matrices
saveRDS(ILV_mat, "ILV_mat.RData")
saveRDS(ILV_mat_int, "ILV_mat_int.RData")


# HRO Matrix ------------------------------------------------------

# Aggregate list into one homerange overlap matrix
list_years <- list_years_int

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
dolph.sp_int <- create_coord_data(list_years_int)

# Use the calculated extent in kernelUD
dolph.sp <- dolph.sp_int
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
saveRDS(kov, "kov_int.RDS")

# GR Matrix ------------------------------------------------------
# HI Matrices ------------------------------------------------------

# Visualize data: HAB v HI
HAB_HI_data <- orig_data[, c("Year", "ConfHI")]
HAB_HI_data$ConfHI <- ifelse(HAB_HI_data$ConfHI != "0", 1, 0)
HAB_HI_data <- aggregate(ConfHI ~ Year, data = HAB_HI_data, FUN = function(x) sum(x == 1))
HAB_HI_data$HAB <- c(rep(0, 4), 5, 0, 12, 8, 18, 5, 38, 19, 2, rep(0, 5))
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
aux_int <- aux_data(list_years_int)

# Categorize ID to Foraging
ID_forg <- function(aux_data) {
  IDbehav <- lapply(aux_data, function(df) {
    df <- table(df$Code, df$Foraging)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df <- df[, c(1, 3)]
    colnames(df) <- c("Code", "Forg_Freq")
    df <- aggregate(. ~ Code, data = df, sum)
    df
  })
  return(IDbehav)
}

IDbehav <- ID_forg(aux)
IDbehav_int <- ID_forg(aux_int)

# Separate HI Behaviors
#' BG = Beg: F, G
#' SD = Scavenge and Depredation: B, C, D, E
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
aux_int <- subset_HI(aux_int)

# Categorize DiffHI to IDs
diff_raw <- function(aux_data) {
  rawHI_diff <- lapply(aux_data, function(df) {
    table_df <- as.data.frame(table(df$Code, df$DiffHI))
    colnames(table_df) <- c("Code", "DiffHI", "Freq")
    return(table_df)
  })}

rawHI_diff <- diff_raw(aux)
rawHI_diff_int <- diff_raw(aux_int)

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

# Two period data
IDbehav_BG <- get_IDHI("BG", IDbehav, rawHI_diff)
IDbehav_FG <- get_IDHI("FG", IDbehav, rawHI_diff)
IDbehav_SD <- get_IDHI("SD", IDbehav, rawHI_diff)
# Three period data
IDbehav_BG_int <- get_IDHI("BG", IDbehav_int, rawHI_diff_int)
IDbehav_FG_int <- get_IDHI("FG", IDbehav_int, rawHI_diff_int)
IDbehav_SD_int <- get_IDHI("SD", IDbehav_int, rawHI_diff_int)

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
rawHI_int <- clump_behav(aux_int)

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
IDbehav_HI_int <- create_IDbehav_HI(IDbehav_int, rawHI_int)

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

# Two period data
prob_HI <- Prop_HI(IDbehav_HI)
prob_BG <- Prop_HI(IDbehav_BG)
prob_SD <- Prop_HI(IDbehav_SD)
prob_FG <- Prop_HI(IDbehav_FG)
# Three period data
prob_HI_int <- Prop_HI(IDbehav_HI_int)
prob_BG_int <- Prop_HI(IDbehav_BG_int)
prob_SD_int <- Prop_HI(IDbehav_SD_int)
prob_FG_int <- Prop_HI(IDbehav_FG_int)

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

# Two period data
dist_HI <- dis_matr(prob_HI, nxn)
dist_BG <- dis_matr(prob_BG, nxn)
dist_SD <- dis_matr(prob_SD, nxn)
dist_FG <- dis_matr(prob_FG, nxn)

saveRDS(dist_HI, "dist_HI.RData")
saveRDS(dist_BG, "dist_BG.RData")
saveRDS(dist_SD, "dist_SD.RData")
saveRDS(dist_FG, "dist_FG.RData")

# Three period data
dist_HI_int <- dis_matr(prob_HI_int, nxn_int)
dist_BG_int <- dis_matr(prob_BG_int, nxn_int)
dist_SD_int <- dis_matr(prob_SD_int, nxn_int)
dist_FG_int <- dis_matr(prob_FG_int, nxn_int)

saveRDS(dist_HI_int, "dist_HI_int.RData")
saveRDS(dist_BG_int, "dist_BG_int.RData")
saveRDS(dist_SD_int, "dist_SD_int.RData")
saveRDS(dist_FG_int, "dist_FG_int.RData")


###########################################################################
# PART 4: Run MCMC GLMM ---------------------------------------------

# Read in social association matrix and listed data
## Two period data
dist_HI <- readRDS("dist_HI.RData") # HI Sim Matrix
dist_BG <- readRDS("dist_BG.RData") # BG Sim Matrix
dist_FG <- readRDS("dist_FG.RData") # FG Sim Matrix
dist_SD <- readRDS("dist_SD.RData") # SD Sim Matrix
ILV_mat <-readRDS("ILV_mat.RData") # Age and Sex Matrices
kov <- readRDS("kov.RDS")  # Home range overlap
nxn <- readRDS("nxn.RData") # Association Matrix
SE_list <- readRDS("SE_list.RData")
## Three period data
dist_HI <- readRDS("dist_HI_int.RData") # HI Sim Matrix
dist_BG <- readRDS("dist_BG_int.RData") # BG Sim Matrix
dist_FG <- readRDS("dist_FG_int.RData") # FG Sim Matrix
dist_SD <- readRDS("dist_SD_int.RData") # SD Sim Matrix
ILV_mat <-readRDS("ILV_mat_int.RData") # Age and Sex Matrices
kov <- readRDS("kov_int.RDS")  # Home range overlap
nxn <- readRDS("nxn_int.RData") # Association Matrix
SE_list <- readRDS("SE_list_int.RData")

# Prepare random effect for MCMC
num_nodes <- lapply(nxn, function(df) dim(df)[1])
node_names <- lapply(nxn, function(df) colnames(df))

# Separate IDs into i and j
node_ids_i <- lapply(num_nodes, function(df) matrix(rep(1:df, each = df), nrow = df, ncol = df))
node_ids_j <- lapply(node_ids_i, function(df) t(df))

# Format data
period = 3
upper_tri <- lapply(nxn, function(df) upper.tri(df, diag = TRUE))
edge_nxn <- abind(lapply(nxn, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
HAB_data <- cbind(c(edge_nxn[,1], edge_nxn[,2]), c(rep(0, nrow(edge_nxn)), rep(1, nrow(edge_nxn)))) # Two
HAB_data <- cbind(c(edge_nxn[,1], edge_nxn[,2], edge_nxn[,3]), c(rep(1, nrow(edge_nxn)), rep(2, nrow(edge_nxn)), rep(3, nrow(edge_nxn)))) # Three
HI <- abind(lapply(dist_HI, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = period)
BG <- abind(lapply(dist_BG, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = period)
FG <- abind(lapply(dist_FG, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = period)
SD <- abind(lapply(dist_SD, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = period)
#SE <- abind(lapply(SE_list, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)
one <- lapply(seq_along(node_ids_i), function(i) factor(as.vector(node_names[[i]][node_ids_i[[i]][upper_tri[[i]]]]), levels = node_names[[i]]))
two <- lapply(seq_along(node_ids_j), function(i) factor(as.vector(node_names[[i]][node_ids_j[[i]][upper_tri[[i]]]]), levels = node_names[[i]]))

# Put data into a dataframe
df_list = data.frame(edge_weight = HAB_data[, 1],
                     HAB = HAB_data[, 2],
                     HRO = unlist(lapply(kov, function (df) df[upper.tri(df, diag = TRUE)])),
                     sex_similarity = rep(ILV_mat[[1]][upper.tri(ILV_mat[[1]], diag = TRUE)], period),
                     age_difference = rep(ILV_mat[[2]][upper.tri(ILV_mat[[2]], diag = TRUE)], period),
                     #GR = gr_list,
                     HI_differences = c(HI[,c(1:period)]),
                     # BG_differences = c(BG[,c(1:period)]),
                     # FG_differences = c(FG[,c(1:period)]),
                     # SD_differences = c(SD[,c(1:period)]),
                     #Obs.Err = c(SE[,1], SE[,2]),
                     node_id_1 = unlist(one),
                     node_id_2 = unlist(two))

# Multimembership models in MCMCglmm
## HI Behavior Combined Two Year Period ##
fit_mcmc.1 <- MCMCglmm(edge_weight ~ HI_differences * HAB + HRO + age_difference + sex_similarity, 
                     random=~mm(node_id_1 + node_id_2), data = df_list, nitt = 20000)
summary(fit_mcmc.1)

fit_mcmc.1 <- MCMCglmm(edge_weight ~ BG_differences * HAB + FG_differences * HAB + SD_differences * HAB + 
                         HRO + age_difference + sex_similarity, 
                       random=~mm(node_id_1 + node_id_2), data = df_list, nitt = 20000)
summary(fit_mcmc.1)

## HI Behavior Combined Three Year Period ##
fit_mcmc.2 <- MCMCglmm(edge_weight ~ HI_differences * HAB + HRO + age_difference + sex_similarity, 
                       random=~mm(node_id_1 + node_id_2), data = df_list, nitt = 20000) 
summary(fit_mcmc.2)

# Check for model convergence
plot(fit_mcmc.2$Sol)
plot(fit_mcmc.2$VCV)

# Extract Posteriors
posterior <- fit_mcmc.2$Sol

# Plot the posterior distribution
mcmc_intervals(posterior, pars = c("(Intercept)", "HI_differences", "HI_differences:HAB", 
                                   "HAB", "age_difference", "sex_similarity", "HRO"))
mcmc_areas(
  posterior, 
  pars = c("(Intercept)", "HI_differences:HAB", "HI_differences", "HAB"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

