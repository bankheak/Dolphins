# Network Analysis of Anthropogenic Influences in Bottlenose dolphins

# Local Network Analysis Hypothesis #2 #

# Set working directory here
setwd("../data")

# Load all necessary packages
library(tnet) # For weights
library(igraph) # Measure centrality here
library(assortnet) # associative indices
library(ggplot2) # Visualization
library(abind) # array
library(MCMCglmm) # MCMC models
library(coda)
library(bayesplot) # plot parameters
library(doParallel)
source("../code/functions.R") # Matrix_to_edge_list

# Read in full datasheet and list (after wrangling steps)
list_years <- readRDS("list_years.RData") # (1995-2000)/(2001-2006)/(2007-2012)
nxn <- readRDS("nxn.RData") # association matrix of list_years

###########################################################################
# PART 1: Wrangle Data ---------------------------------------------

# Read in full datasheet and list
list_years <- readRDS("list_years.RData") # (1995-2000)/(2001-2006)/(2007-20012)
nxn <- readRDS("nxn.RData") # association matrix of list_years

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

# Look at how many individuals have HI
HI_1 <- unique(aux[[1]]$Code[aux[[1]]$DiffHI != "None"])
HI_2 <- unique(aux[[2]]$Code[aux[[2]]$DiffHI != "None"])
HI_3 <- unique(aux[[3]]$Code[aux[[3]]$DiffHI != "None"])
length(unique(c(HI_1, HI_2, HI_3)))

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

IDbehav_BG <- get_IDHI("BG", IDbehav, rawHI_diff)
IDbehav_FG <- get_IDHI("FG", IDbehav, rawHI_diff)
IDbehav_SD <- get_IDHI("SD", IDbehav, rawHI_diff)

# Get total number of HI individuals
BG_IDs <- unique(unlist(lapply(IDbehav_BG, function (df) unique(df$Code[df$HI > 0]))))
FG_IDs <- unique(unlist(lapply(IDbehav_FG, function (df) unique(df$Code[df$HI > 0]))))
SD_IDs <- unique(unlist(lapply(IDbehav_SD, function (df) unique(df$Code[df$HI > 0]))))
ovrlap_IDs <- intersect(intersect(BG_IDs, FG_IDs), SD_IDs)


###########################################################################
# PART 2: Calculate Local Metrics ---------------------------------------------

# Edgelist: Nodes (i & j) and edge (or link) weight
n.cores <- detectCores()
registerDoParallel(n.cores)
el_years <- lapply(nxn, function (list) matrix_to_edgelist(list, rawdata = FALSE, idnodes = FALSE))
### End parallel processing
stopImplicitCluster()

# Save the el_list
saveRDS(el_years, "el_years.RData")
el_years <- readRDS("el_years.RData")

# Set the node names based on row names
get_names <- function (matrix, metric) {
  row_names <- lapply(matrix, function (df) {rownames(df)})
  for (i in seq_along(metric)) {
    metric[[i]][,1] <- row_names[[i]]
  }
  return(metric)
}

# Betweenness centrality 
between <- lapply(el_years, function (df) {betweenness_w(df, alpha = 1)})
between_diffs <- get_names(nxn, between)

compare_between <- merge(
  merge(between_diffs[[1]], between_diffs[[2]], by = "node"),
  between_diffs[[3]], by = "node"
)

colnames(compare_between) <- c("ID", "Before_HAB", "During_HAB", "After_HAB")
compare_between[, c(2:4)] <- sapply(compare_between[, c(2:4)], as.numeric)

# Degree and strength centrality 
strength <- lapply(el_years, function (df) {degree_w(df, measure=c("degree","output"), type="out", alpha=1)})
strength_diffs <- get_names(nxn, strength)

compare_strength <- merge(
  merge(strength_diffs[[1]], strength_diffs[[2]], by = "node"),
  strength_diffs[[3]], by = "node"
)

colnames(compare_strength) <- c("ID", "Before_HAB_degree", "Before_HAB_strength", 
                                "During_HAB_degree", "During_HAB_strength", 
                                "After_HAB_degree", "After_HAB_strength")

compare_strength[, c(2:7)] <- sapply(compare_strength[, c(2:7)], as.numeric)

# Look at all of the local metrics together
HI_data <-  subset_HI(list_years)

## Add a column containing HI type
names_BG <- lapply(HI_data, function (df) {
  as.vector(unique(df$Code[df$DiffHI == "BG"]))})
names_SD <- lapply(HI_data, function (df) {
  as.vector(unique(df$Code[df$DiffHI == "SD"]))})
names_FG <- lapply(HI_data, function (df) {
  as.vector(unique(df$Code[df$DiffHI == "FG"]))})

### Find all of the IDs that match to make sure below works
matching_unique_ids <- list()
for (i in 1:3) {
  matching_unique_ids[[i]] <- unique(c(names_BG[[i]], names_FG[[i]], names_SD[[i]]))
}

names_NF <- list()
for (i in 1:3) {
  unique_codes <- unique(HI_data[[i]]$Code)
  
# Check for codes that are not in any of the names_BG, names_FG, names_SD
names_NF[[i]] <- unique_codes[!(unique_codes %in% names_BG[[i]] | 
                                    unique_codes %in% names_FG[[i]] | 
                                    unique_codes %in% names_SD[[i]])]}


HI_list <- list(BG = names_BG, FG = names_FG, SD = names_SD, NF = names_NF)

# Combine the data
local_metrics_HI <- data.frame(ID = compare_between$ID,
                               Period = c(rep("1-Before_HAB", nrow(compare_between)), 
                                          rep("2-During_HAB", nrow(compare_between)),
                                          rep("3-After_HAB", nrow(compare_between))),
                               Between = c(compare_between$Before_HAB, 
                                           compare_between$During_HAB,
                                           compare_between$After_HAB),
                               Degree = c(compare_strength$Before_HAB_degree, 
                                          compare_strength$During_HAB_degree,
                                          compare_strength$After_HAB_degree),
                               Strength = c(compare_strength$Before_HAB_strength, 
                                            compare_strength$During_HAB_strength,
                                            compare_strength$After_HAB_strength))

# Add HI_type column
local_metrics_HI_1 <- local_metrics_HI[local_metrics_HI$Period == "1-Before_HAB", ]
local_metrics_HI_2 <- local_metrics_HI[local_metrics_HI$Period == "2-During_HAB", ]
local_metrics_HI_3 <- local_metrics_HI[local_metrics_HI$Period == "3-After_HAB", ]
local_list <- list(local_metrics_HI_1, local_metrics_HI_2, local_metrics_HI_3)

## Initialize a new dataframe to store the results
result_df <- data.frame()
## Initialize a counter
for (p in 1:3) {
  counter <- 0
  result_df_new <- data.frame()
  for (i in HI_list) {
    # Increment the counter
    counter <- counter + 1
    index <- local_list[[p]]$ID %in% i[[p]]
    # Create a new row for each ID that falls into four categories
    new_rows <- local_list[[p]][index, ]
    new_rows$HI <- names(HI_list[counter])
    # Append the new rows to the result dataframe
    result_df_new <- rbind(result_df_new, new_rows)
  }
  result_df <- rbind(result_df, result_df_new)
}

# Make period a binary variable and HI a factorial
result_df$During <- ifelse(result_df$Period == "2-During_HAB", 1, 0)
result_df$After <- ifelse(result_df$Period == "3-After_HAB", 1, 0)
result_df$HI <- as.factor(result_df$HI)

# Save dataset
saveRDS(result_df, "result_df.RData")


###########################################################################
# PART 3: Run Model ---------------------------------------------

# Set up data
result_df <- readRDS("result_df.RData")

## Between
ggplot(result_df, aes(x = Period, y = Between, fill = HI)) + 
  geom_boxplot()
## Strength
ggplot(result_df, aes(x = Period, y = Strength, fill = HI)) + 
  geom_boxplot() 
## Degree
ggplot(result_df, aes(x = Period, y = Degree, fill = HI)) + 
  geom_boxplot() 
  
# Check distributions
hist(result_df$Between) # continuous
hist(result_df$Degree)
hist(result_df$Strength)

# Make dummy variables
result_df$BG <- ifelse(result_df$HI == "BG", 1, 0)
result_df$FG <- ifelse(result_df$HI == "FG", 1, 0)
result_df$SD <- ifelse(result_df$HI == "SD", 1, 0)

# CHeck if just HI behavior changes over time
# Try and plot multilayer network
# Plot the x and y and plot netwirks in columns and rows
# In the same plot graph dist of centrality metrics

# Look into nodal regression
## HI Behavior Combined Two Year Period ##
fit_mcmc.b <- MCMCglmm(Between ~ BG * During + FG * During + SD * During +
                       BG * After + FG * After + SD * After, data = result_df, nitt = 10000, family = "poisson")
summary(fit_mcmc.b) # Might need to use a negative binomial dist

fit_mcmc.s <- MCMCglmm(Strength ~ BG * During + FG * During + SD * During +
                       BG * After + FG * After + SD * After, data = result_df, nitt = 10000)
summary(fit_mcmc.s)

fit_mcmc.d <- MCMCglmm(Degree ~ BG * During + FG * During + SD * During +
                       BG * After + FG * After + SD * After, data = result_df, nitt = 10000)
summary(fit_mcmc.d)

# Check for model convergence
model <- fit_mcmc.d
plot(model$Sol)
plot(model$VCV)

# Extract Posteriors
posterior <- model$Sol

# Plot the posterior distribution
mcmc_intervals(posterior, pars = c("(Intercept)", "BG", "FG", "SD",
                                   "During", "After", "BG:During", "During:FG", "During:SD",
                                   "BG:After", "FG:After", "SD:After"))
mcmc_areas(
  posterior, 
  pars = c("(Intercept)", "BG", "FG", "SD",
           "During", "After", "BG:During", "During:FG", "During:SD",
           "BG:After", "FG:After", "SD:After"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

# Test if model is good for predicting data
# Make empty list for each row's distribution
edge_weight <- Obs.edge_weight <- vector("list", length = nrow(df_list))
# Make an empty vector for the true and false values
exp.obs <- NULL
posterior <- as.data.frame(posterior)

for (i in 1:nrow(df_list)) {
  # Expected bill length
  edge_weight[[i]] <-  posterior[,"(Intercept)"] + posterior[,"HI_differences:HAB"]*(df_list$HI_differences[i] * df_list$HAB[i]) + 
    posterior[,"HI_differences"] * df_list$HI_differences[i] + posterior[,"HAB"] * df_list$HAB[i]
  posterior[,"HRO"]*df_list$HRO[i] + 
    posterior[,"age_difference"]*df_list$age_difference[i] + 
    posterior[,"sex_similarity"]*df_list$sex_similarity[i]
  
  # Observed bill length
  Obs.edge_weight[[i]] <- rnorm(n = 1700, mean = edge_weight[[i]], sd = rep(sd(edge_weight[[i]]), nrow(df_list)))
  
  # Calculate how often observed values fall into expected
  exp.obs[i] <- df_list$edge_weight[i] >= quantile(Obs.edge_weight[[i]], c(0.025, 0.975))[1] & 
    df_list$edge_weight[i] <= quantile(Obs.edge_weight[[i]], c(0.025, 0.975))[2]
}

sum(exp.obs)/length(exp.obs)

###########################################################################
# PART 4: Circular heat map ---------------------------------------------

# Set up data
result_df <- readRDS("result_df.RData")

none_df <- result_df[!(result_df$Period == 1 & result_df$HI == "NF"), ]
result_df <- none_df[!(none_df$Period == 0 & !(none_df$ID %in% none_df$ID[none_df$Period == 1]) & none_df$HI == "NF"), ]
ID_1 <- unique(result_df$ID[result_df$Period==0])
ID_2 <- unique(result_df$ID[result_df$Period==1])
ID_length <- rep(c(ID_1, ID_2), 3)
Period <- c("Pre-HAB", "Post_HAB")

# Filter data for Between and unique IDs
unique_data_b <- result_df[!duplicated(result_df[, c('ID', 'Period')]), ]

# Combine the sets of data
B_data <- data.frame(ID = ID_length, 
                     Period = rep(c(rep(Period[1], length(ID_1)), 
                                rep(Period[2], length(ID_2))), 3),
                     HI = rep(unique_data_b$HI, 3),
                     Value = c(c(scale(c(unique_data_b$Between[unique_data_b$Period == 0], 
                                            unique_data_b$Between[unique_data_b$Period == 1]))),
                             c(scale(c(unique_data_b$Strength[unique_data_b$Period == 0], 
                                          unique_data_b$Strength[unique_data_b$Period == 1]))),
                             c(scale(c(unique_data_b$Degree[unique_data_b$Period == 0], 
                                        unique_data_b$Degree[unique_data_b$Period == 1])))),
                     Metric = c(rep("Betweeness", length(ID_length)),
                                rep("Strength", length(ID_length)),
                                rep("Degree", length(ID_length))))

# plotting the heatmap
heatmap_list <- list()
count <- 1
for (j in c("BG", "FG", "SD")) {
    HI_data <- B_data[B_data$HI == j, ]
    matched_ids <- B_data[B_data$ID %in% 
                            intersect(B_data$ID[B_data$Period == "Pre-HAB" & B_data$HI == "NF"], 
                                      B_data$ID[B_data$Period == "Post_HAB" & B_data$HI == j]), ]
    matched_ids <- matched_ids[matched_ids$HI == "NF", ]
    subset_data <- rbind(HI_data, matched_ids)
    
    heatmap_list[[count]] <- ggplot(subset_data, aes(ID, Period, fill = Value)) + 
      geom_tile() + 
      theme_minimal() + 
      scale_fill_gradient(low="white", high="red") + 
      labs(title = "Heatmap of Centrality Metrics of Dolphins", 
           x ="IDs", y ="Period")
    count <- 1 + 1
  }

# Plot heatmaps
heatmap_list[[1]] # BG
heatmap_list[[2]] # FG
heatmap_list[[3]] # SD
