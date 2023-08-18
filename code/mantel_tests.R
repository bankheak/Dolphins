# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# MRQAP TESTS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: Create HI Disimilarity Matrix  ------------------------------------------------

## load all necessary packages
library(ade4) # Look at Dai Shizuka/Jordi Bascompte
library(ncf) # For weights
library(vegan)
library(igraph) # graph_adj
require(asnipe) # mrqap.dsp

# Read file in to retain ILV
sample_data <- read.csv("sample_data.csv")
kov <- readRDS("kov.RDS")  # Home range overlap
kov <- as.dist(kov)

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

# HI behaviors should be partitioned into 3 different types
#' B = Begging (direct provisioning): F, G, H
#' P = patrolling/scavenging (indirect): A, B, C
#' D = foraging around fixed gear (humans not present):D, E, P
# Fix the code using ifelse statements
for (i in seq_along(aux)) {
  
  aux[[i]]$ConfHI <- ifelse(aux[[i]]$ConfHI %in% c("F", "G", "H"), "B",
                            ifelse(aux[[i]]$ConfHI %in% c("A", "B", "C", "D", "E"), "S", 
                                   ifelse(aux[[i]]$ConfHI %in% c("P"), "D", "0")))
  
}

# Categorize ConfHI to IDs
rawHI <- lapply(aux, function(df) {
  df <- as.matrix(table(df$Code, df$ConfHI))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  colnames(df) <- c("Code", "ConfHI", "Freq")
  df
})

# Create a different frequency count for each HI behavior
get_IDHI <- function(confHI) {
  lapply(seq_along(IDbehav), function(i) {
    df <- IDbehav[[i]]
    df$HI <- rawHI[[i]]$Freq[rawHI[[i]]$ConfHI == confHI & rawHI[[i]]$ConfHI != "0"]
    colnames(df) <- c("Code", "Foraging", "HI")
    df
  })
}

IDbehav_Beg <- get_IDHI("B")
IDbehav_Pat <- get_IDHI("S")
IDbehav_Dep <- get_IDHI("D")

saveRDS(IDbehav_Beg, file = "../data/IDbehav_Beg.RData")
saveRDS(IDbehav_Pat, file = "../data/IDbehav_Pat.RData")
saveRDS(IDbehav_Dep, file = "../data/IDbehav_Dep.RData")

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

prob_Beg <- Prop_HI(IDbehav_Beg)
prob_Pat <- Prop_HI(IDbehav_Pat)
prob_Dep <- Prop_HI(IDbehav_Dep)

# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
dis_matr <- function(IDbehav) {
  dissimilarity_HI <- list()
  for (i in seq_along(IDbehav)) {
    fake_HIprop <- IDbehav[[i]]$HIprop
    dissimilarity_HI[[i]] <- as.matrix(dist(matrix(fake_HIprop), method = "euclidean"))
    dissimilarity_HI[[i]][is.na(dissimilarity_HI[[i]])] <- 0
    dissimilarity_HI[[i]] <- as.dist(dissimilarity_HI[[i]]) # HI dissimilarity
  }
  dissimilarity_HI
}

dist_Beg <- dis_matr(prob_Beg)
dist_Pat <- dis_matr(prob_Pat)
dist_Dep <- dis_matr(prob_Dep)

###########################################################################
# PART 2: Run Mantel Tests  ------------------------------------------------

# Dissimilarity Mantel Test
year <- 5
HI_test <- mantel.rtest(dolp_dist[[year]], dist_Pat[[year]], nrepet = 1000)
plot(HI_test)
# So far no correlation with HI engagement and associations

## HRO
hro_test <- mantel.rtest(dolp_dist, kov, nrepet = 1000)
plot(hro_test)


###########################################################################
# PART 3: Create MRQAP Models  ------------------------------------------------

# Set a number of permutations
Nperm <- 1000

# Calculate QAP correlations for the association response matrix
mrqap <- mrqap.dsp(nxn[[year]] ~ dist_Beg[[year]] + kov + sex + age,
                   randomisations = Nperm,
                   intercept = FALSE,
                   test.statistic = "beta")

