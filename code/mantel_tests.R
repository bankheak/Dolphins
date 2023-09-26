# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# MRQAP TESTS
###########################################################################

# Set working directory here
setwd("../data")

###########################################################################
# PART 1: Create HI Disimilarity Matrix  ------------------------------------------------

## load all necessary packages
library(ade4) # Look at Dai Shizuka/Jordi Bascompte
require(asnipe) # mrqap.dsp
library(assortnet) # associative indices
library(ggplot2) # Visualization
library(doParallel) # For faster coding
library(MCMCglmm) # MCMC models
library(brms) # Baysian

# Read in social association matrix and listed data
kov <- readRDS("kov.RDS")  # Home range overlap
nxn <- readRDS("nxn.RData") # Association Matrix
list_years <- readRDS("list_years.RData") # Data listed into periods

# Read file in to retain ILV
kov_sexage <- readRDS("kov_sexage.RDS")  
nxn_sexage <- readRDS("nxn_sexage.RData")
list_sexage_years <- readRDS("list_sexage_years.RData")

# Read file in to retain only HI IDs
kov_HI <- readRDS("kov_HI.RDS") 
nxn_HI <- readRDS("nxn_HI.RData") 

# Transforming SRI similarity into distance
dolp_dist <- lapply(nxn, function(df) {
  df + 0.00001
  1 - df
  ## Remove the redundant cells and the diagonal 
  as.dist(df)
})

# Sex similarity matrix
sex_list <- lapply(list_sexage_years, function(df) {
  
  ## Empty matrix to store sex similarity
  num_ID <- length(unique(df$Code))
  sex_matrix <- matrix(NA, nrow = num_ID, ncol = num_ID, dimnames = list(unique(df$Code), unique(df$Code)))
  
  # Fill in similarity of sex
  for (i in 1:num_ID) {
    for (j in 1:num_ID) {
      if (df$Sex[i] == df$Sex[j]) {
        sex_matrix[i, j] <- 1  # Same sex
      } else {
        sex_matrix[i, j] <- 0  # Different sex
      }
    }
  }
  
  return(sex_matrix)
})

# Age similarity matrix
age_list <- lapply(list_sexage_years, function(df) {
  
  ## Empty matrix to store sex similarity
  num_ID <- length(unique(df$Code))
  age_matrix <- matrix(NA, nrow = num_ID, ncol = num_ID, dimnames = list(unique(df$Code), unique(df$Code)))
  
  # Fill in similarity of sex
  for (i in 1:num_ID) {
    for (j in 1:num_ID) {
      age_matrix[i, j] <- abs(df$Age[i] - df$Age[j])
    }
  }
  return(age_matrix)
})

# Extract specific columns from each data frame in list_years
aux_data <- function(list_years) {
aux <- lapply(list_years, function(df) {
  data.frame(
    Code = df$Code,
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
aux_sexage <- aux_data(list_sexage_years)
aux_HI <- aux_data(list_HI_years)

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
IDbehav_sexage <- ID_forg(aux_sexage)
IDbehav_HI <- ID_forg(aux_HI)

# HI behaviors should be partitioned into 3 different types---------------------
#' B = Beg: F, G
#' P = Scavenge and Depredation: B, C, D, E
#' D = Fixed Gear Interaction: P
# Change the code using ifelse statements
subset_HI <- function(aux_data) {
  for (i in seq_along(aux_data)) {
    aux_data[[i]]$DiffHI <- ifelse(aux_data[[i]]$ConfHI %in% c("F", "G"), "BG",
                                   ifelse(aux_data[[i]]$ConfHI %in% c("B", "C", "D", "E"), "SD",
                                          ifelse(aux_data[[i]]$ConfHI %in% c("P"), "FG", "None")))
  }
  return(aux_data)  # Return the modified list of data frames
}

aux <- subset_HI(aux)
aux_HI <- subset_HI(aux_HI)

# Categorize DiffHI to IDs
diff_raw <- function(aux_data) {
rawHI_diff <- lapply(aux_data, function(df) {
  table_df <- as.data.frame(table(df$Code, df$DiffHI))
  colnames(table_df) <- c("Code", "DiffHI", "Freq")
  return(table_df)
})}

rawHI_diff <- diff_raw(aux)
rawHI_diff_HI <- diff_raw(aux_HI)

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
IDbehav_SD <- get_IDHI("SD", IDbehav, rawHI_diff)
IDbehav_FG <- get_IDHI("FG", IDbehav, rawHI_diff)

# Not including zeros
IDbehav_BG <- get_IDHI("BG", IDbehav_HI, rawHI_diff_HI)
IDbehav_SD <- get_IDHI("SD", IDbehav_HI, rawHI_diff_HI)
IDbehav_FG <- get_IDHI("FG", IDbehav_HI, rawHI_diff_HI)

# Clump all the HI behaviors together------------------------------------------
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
rawHI_sexage <- clump_behav(aux_sexage)
rawHI_HI <- clump_behav(aux_HI)

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
IDbehav_HI_sexage <- create_IDbehav_HI(IDbehav_sexage, rawHI_sexage)
IDbehav_HI_HI <- create_IDbehav_HI(IDbehav_HI, rawHI_HI)

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
prop_HI_sexage <- Prop_HI(IDbehav_HI_sexage)
prop_HI_HI <- Prop_HI(IDbehav_HI_HI)
prob_Beg <- Prop_HI(IDbehav_Beg)
prob_Pat <- Prop_HI(IDbehav_Pat)
prob_Dep <- Prop_HI(IDbehav_Dep)

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
dist_HI_sexage <- dis_matr(prop_HI_sexage)
dist_HI_HI <- dis_matr(prop_HI_HI)
dist_Beg <- dis_matr(prob_Beg)
dist_Pat <- dis_matr(prob_Pat)
dist_Dep <- dis_matr(prob_Dep)


###########################################################################
# PART 2: Create MRQAP Models  ------------------------------------------------

# Check for collinearity 
# Check if it is based off zeros
## Create a list of predictor matrices
predictor_matrices <- list(dist_HI[[year]], dist_Beg[[year]], 
                           dist_Pat[[year]], dist_Dep[[year]])

## Calculate correlation matrix
num_predictors <- length(predictor_matrices)
correlation_matrix <- matrix(NA, nrow = num_predictors, ncol = num_predictors)

for (i in 1:num_predictors) {
  for (j in 1:num_predictors) {
    mtest <- mantel.rtest(as.dist(predictor_matrices[[i]]), as.dist(predictor_matrices[[j]]), nrepet=999)
    correlation_matrix[i, j] <- mtest$obs
  }
}

## Print the correlation matrix
print(correlation_matrix) # It seems that BEG and HI are highly correlated

# Set a number of permutations and year
year <- 5
Nperm <- 1000

# Calculate QAP correlations for the association response matrix

## Without sex and age included
mrqap_full <- mrqap.dsp(nxn[[year]] ~ kov[[year]] + dist_HI[[year]],
                   randomisations = Nperm,
                   intercept = FALSE,
                   test.statistic = "beta")

## Without sex and age included and with behaviors divided
mrqap_sepHI <- mrqap.dsp(nxn[[year]] ~ kov[[year]] +
                     dist_Beg[[year]] + dist_Dep[[year]] + dist_Pat[[year]],
                   randomisations = Nperm,
                   intercept = FALSE,
                   test.statistic = "beta")

## With sex and age included
mrqap_sexage <- mrqap.dsp(nxn_sexage[[year]] ~ kov_sexage[[year]] + 
                            sex_list[[year]] + age_list[[year]] + 
                            dist_HI_sexage[[year]],
                   randomisations = Nperm,
                   intercept = FALSE,
                   test.statistic = "beta")

## With only HI individuals included
mrqap_HIonly <- mrqap.dsp(nxn_HI[[year]] ~ kov_HI + dist_HI_HI[[year]],
                        randomisations = Nperm,
                        intercept = FALSE,
                        test.statistic = "beta")


###########################################################################
# PART 3: Create MCMC GLMMs  ------------------------------------------------

# Prepare dataframe for MCMC
num_nodes <- lapply(nxn_sexage, function(df) {dim(df)[1]})

## Seperate IDs into i and j
node_ids_i <- lapply(num_nodes, function(df) {matrix(rep(1:df, df), df, df)})
node_ids_j <- lapply(node_ids_i, function(df) {t(df)})

df_list <- list()  
for (i in seq_along(nxn_sexage)) {
  df_dolp <- data.frame(edge_weight = nxn_sexage[[i]][upper.tri(nxn_sexage[[i]])], 
                        age_difference = age_list[[i]][upper.tri(age_list[[i]])],
                        sex_difference = sex_list[[i]][upper.tri(sex_list[[i]])],
                        HI_differences = dist_HI_sexage[[i]][upper.tri(dist_HI_sexage[[i]])],
                        HRO = kov_sexage[[i]][upper.tri(kov_sexage[[i]])],
                        node_id_1 = factor(node_ids_i[[i]][upper.tri(node_ids_i[[i]])],
                                           levels = 1:num_nodes[[i]]),
                        node_id_2 = factor(node_ids_j[[i]][upper.tri(node_ids_j[[i]])], 
                                           levels = 1:num_nodes[[i]]))
  df_list[[i]] <- df_dolp
}

# Multimembership models in MCMCglmm
year <- 5
fit_mcmc <- MCMCglmm(edge_weight ~ HI_differences + HRO + age_difference + sex_difference, 
                     random=~mm(node_id_1 + node_id_2), data=df_list[[year]])
summary(fit_mcmc)


###########################################################################
# PART 4: Assortivity Index Based on HI Over Time  ------------------------------------------------

# Match Code with matrix and vector
get_HI_vector <- function(prop_HI) {
  HI_vector <- lapply(seq_along(nxn), function(i) {
    matrix_index <- match(rownames(nxn[[i]]), prop_HI[[i]]$Code)
    reordered_prob_HI <- prop_HI[[i]][matrix_index, ]
    return(reordered_prob_HI)
  })
  return(HI_vector)
}

# Get each combined and seperate HI
HI_vector <- get_HI_vector(prob_HI)
Beg_vector <- get_HI_vector(prob_Beg)
Pat_vector <- get_HI_vector(prob_Pat)
Dep_vector <- get_HI_vector(prob_Dep)

# Look at HI assortivity coefficient over periods
calculate_assortment <- function(HI_vector) {
  n.cores <- detectCores()
  registerDoParallel(n.cores)
  
  assort_HI <- NULL
  # se <- NULL
  for (i in seq_along(nxn)) {
    coeff <- assortment.continuous(nxn[[i]], HI_vector[[i]][, "HIprop"], SE = FALSE)
    assort_HI[i] <- coeff$r
    # se[i] <- coeff$se
  }
  
  # End parallel processing
  stopImplicitCluster()
  
  assort_HI_df <- data.frame(HI_assort = unlist(assort_HI), Year = c(1:7))
  return(assort_HI_df)
}

# Look at HI combined and separate
assort_HI <- calculate_assortment(HI_vector)
assort_Beg <- calculate_assortment(Beg_vector)
assort_Pat <- calculate_assortment(Pat_vector)
assort_Dep <- calculate_assortment(Dep_vector)

# Combine the assort dataframes and add a behavior column
assort_Beg$Behavior <- "Beg"
assort_Pat$Behavior <- "Pat"
assort_Dep$Behavior <- "Dep"

combined_assort <- rbind(assort_Beg, assort_Pat, assort_Dep)
combined_assort$HI_assort <- ifelse(is.na(combined_assort$HI_assort), 0.75, combined_assort$HI_assort)

# Create the combined plot with facets
ggplot(combined_assort, aes(x = Year, y = HI_assort)) +
  geom_point() +
  geom_line() +
  labs(x = "Period", y = "HI assortment") +
  ggtitle("Whisker Plot of HI assortment") +
  theme_minimal() +
  facet_grid(Behavior ~ ., scales = "free_y", space = "free_y")
