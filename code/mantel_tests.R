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
library(ncf) # For weights
library(vegan)
library(igraph) # graph_adj
require(asnipe) # mrqap.dsp
library(assortnet)
library(ggplot2)
library(doParallel)

# Read file in to retain ILV
sample_data <- read.csv("sample_data.csv")
kov <- readRDS("kov.RDS")  # Home range overlap

# Read in social association matrix and data
nxn <- readRDS("nxn.RData")
list_years <- readRDS("list_years.RData")

# Transforming SRI similarity into distance
dolp_dist <- lapply(nxn, function(df) {
  df + 0.00001
  1 - df
  ## Remove the redundant cells and the diagonal 
  as.dist(df)
})

# Sex similarity matrix
sex_list <- lapply(list_years, function(df) {
  
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
age_list <- lapply(list_years, function(df) {
  
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
aux <- lapply(list_years, function(df) {
  data.frame(
    Code = df$Code,
    Behaviors = df$Behaviors,
    HumanInteraction = df$HumanInteraction,
    ConfHI = df$ConfHI
  )
})

# Add the 'Foraging' variable to each data frame in the 'aux' list
aux <- lapply(aux, function(df) {
  df$Foraging <- "Other"
  df$Foraging[grepl(pattern = 'Feed', x = df$Behaviors, ignore.case = FALSE)] <- "Feed"
  df
})

# Categorize ID to Foraging
IDbehav <- lapply(aux, function(df) {
  df <- table(df$Code, df$Foraging)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- df[, c(1, 3)]
  colnames(df) <- c("Code", "Forg_Freq")
  df <- aggregate(. ~ Code, data = df, sum)
  df
})

# HI behaviors should be partitioned into 3 different types---------------------
#' B = Beg: F, G, H
#' P = Patrol: A, B, C
#' D = Depredation: D, E, P
# Change the code using ifelse statements
for (i in seq_along(aux)) {

  aux[[i]]$DiffHI <- ifelse(aux[[i]]$ConfHI %in% c("F", "G", "H"), "Beg",
                            ifelse(aux[[i]]$ConfHI %in% c("A", "B", "C"), "Pat",
                                   ifelse(aux[[i]]$ConfHI %in% c("P", "D", "E"), "Dep", "0")))

}

# Categorize DiffHI to IDs
rawHI_diff <- lapply(aux, function(df) {
  table_df <- as.data.frame(table(df$Code, df$DiffHI))
  colnames(table_df) <- c("Code", "DiffHI", "Freq")
  return(table_df)
})

# Create a frequency count for each HI behavior
get_IDHI <- function(HI) {
  lapply(seq_along(IDbehav), function(i) {
    df <- IDbehav[[i]]
    HI_freq <- rawHI_diff[[i]]$Freq[rawHI_diff[[i]]$DiffHI == HI]
    df$HI <- HI_freq[match(df$Code, rawHI_diff[[i]]$Code)]
    colnames(df) <- c("Code", "Foraging", "HI")
    df
  })
}

IDbehav_Beg <- get_IDHI("Beg")
IDbehav_Pat <- get_IDHI("Pat")
IDbehav_Dep <- get_IDHI("Dep")

saveRDS(IDbehav_Beg, file = "../data/IDbehav_Beg.RData")
saveRDS(IDbehav_Pat, file = "../data/IDbehav_Pat.RData")
saveRDS(IDbehav_Dep, file = "../data/IDbehav_Dep.RData")

# Clump all the HI behaviors together------------------------------------------
for (i in seq_along(aux)) {
aux[[i]]$ConfHI <- ifelse(aux[[i]]$ConfHI != "0", 1, 0)}

# Categorize ConfHI to IDs
rawHI <- lapply(aux, function(df) {
  # Sum up the frequencies of HI by code
  aggregated_df <- aggregate(ConfHI ~ Code, data = df, sum)
  unique_codes_df <- data.frame(Code = unique(df$Code))
  # Merge the unique codes data frame with the aggregated data frame
  merged_df <- merge(unique_codes_df, aggregated_df, by = "Code", all.x = TRUE)
  # Fill missing Freq values (if any) with 0
  merged_df$ConfHI[is.na(merged_df$ConfHI)] <- 0
  return(merged_df)
})

# Get HI Freq
IDbehav_HI <- lapply(seq_along(IDbehav), function(i) {
      df <- IDbehav[[i]]
      df$HI <- rawHI[[i]]$ConfHI
      colnames(df) <- c("Code", "Foraging", "HI")
      df
    })

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
prob_Beg <- Prop_HI(IDbehav_Beg)
prob_Pat <- Prop_HI(IDbehav_Pat)
prob_Dep <- Prop_HI(IDbehav_Dep)

# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
dis_matr <- function(Prop_HI) {
  dissimilarity_HI <- list()
  for (i in seq_along(Prop_HI)) {
    fake_HIprop <- Prop_HI[[i]]$HIprop
    dissimilarity_HI[[i]] <- as.matrix(dist(matrix(fake_HIprop), method = "euclidean"))
    dissimilarity_HI[[i]][is.na(dissimilarity_HI[[i]])] <- 0
    #dissimilarity_HI[[i]] <- as.dist(dissimilarity_HI[[i]]) # HI dissimilarity
  }
  dissimilarity_HI
}

dist_HI <- dis_matr(prob_HI)
dist_Beg <- dis_matr(prob_Beg)
dist_Pat <- dis_matr(prob_Pat)
dist_Dep <- dis_matr(prob_Dep)


###########################################################################
# PART 2: Create MRQAP Models  ------------------------------------------------

# Set a number of permutations and year
year <- 5
Nperm <- 1000

# Calculate QAP correlations for the association response matrix
mrqap <- mrqap.dsp(nxn[[year]] ~ kov[[year]] + dist_HI[[year]] +
                     dist_Beg[[year]] + dist_Dep[[year]] + dist_Pat[[year]],
                   randomisations = Nperm,
                   intercept = FALSE,
                   test.statistic = "beta")

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
