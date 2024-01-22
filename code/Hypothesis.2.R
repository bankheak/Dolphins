# Network Analysis of Anthropogenic Influences in Bottlenose dolphins

# Local Network Analysis Hypothesis #2 #

# Set working directory here
setwd("../data")

# Load all necessary packages
library(tnet) # For weights
library(igraph) # Measure centrality here
library(assortnet) # associative indices
library(kinship2) # genetic relatedness
library(ggplot2) # Visualization
library(abind) # array
library(MCMCglmm) # MCMC models
library(coda)
library(bayesplot) # plot parameters
library(doParallel)
source("../code/functions.R") # Matrix_to_edge_list

###########################################################################
# PART 1: Wrangle Data ---------------------------------------------

# Read in full datasheet and list
list_years <- readRDS("list_years.RData") # (1998-2004)/(2005-2014)
list_years_int <- readRDS("list_years_int.RData") # (1995-2000)/(2001-2006)/(2007-20012)
nxn <- readRDS("nxn.RData") # association matrix of list_years
nxn_int <- readRDS("nxn_int.RData") # association matrix of list_years_int

# Look at how many individuals have HI
HI_1 <- unique(list_years[[1]]$Code[list_years[[1]]$ConfHI != 0])
HI_2 <- unique(list_years[[2]]$Code[list_years[[2]]$ConfHI != 0])
length(HI_1) + length(HI_2)


###########################################################################
# PART 2: Calculate Local Metrics ---------------------------------------------

# Edgelist: Nodes (i & j) and edge (or link) weight
n.cores <- detectCores()
registerDoParallel(n.cores)
el_years <- lapply(nxn, function (list) matrix_to_edgelist(list, rawdata = FALSE, idnodes = FALSE))
### End parallel processing
stopImplicitCluster()

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
  between_diffs[[1]], 
  between_diffs[[2]], 
  by.x = "node", 
  by.y = "node"
)

colnames(compare_between) <- c("ID", "Before.HAB", "After.HAB")
compare_between[, c(2, 3)] <- sapply(compare_between[, c(2, 3)], as.numeric)

# Degree and strength centrality
strength <- lapply(el_years, function (df) {degree_w(df, measure=c("degree","output"), type="out", alpha=1)})
strength_diffs <- get_names(nxn, strength)

compare_strength <- merge(
  strength_diffs[[1]], 
  strength_diffs[[2]], 
  by.x = "node", 
  by.y = "node"
)

colnames(compare_strength) <- c("ID", "Before.HAB_degree", "Before.HAB_strength", "After.HAB_degree", "After.HAB_strength")
compare_strength[, c(2:5)] <- sapply(compare_strength[, c(2:5)], as.numeric)

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
for (i in 1:2) {
  matching_unique_ids[[i]] <- unique(c(names_BG[[i]], names_FG[[i]], names_SD[[i]]))
}

names_NF <- list()
for (i in 1:2) {
  unique_codes <- unique(HI_data[[i]]$Code)
  
  # Check for codes that are not in any of the names_BG, names_FG, names_SD
  names_NF[[i]] <- unique_codes[!(unique_codes %in% names_BG[[i]] | 
                                    unique_codes %in% names_FG[[i]] | 
                                    unique_codes %in% names_SD[[i]])]}


HI_list <- list(BG = names_BG, FG = names_FG, SD = names_SD, NF = names_NF)

# Combine the data
local_metrics_HI <- data.frame(ID = compare_between$ID,
                               Period = c(rep("1-Before.HAB", nrow(compare_between)), rep("2-After.HAB", nrow(compare_between))),
                               Between = c(compare_between$Before.HAB, compare_between$After.HAB),
                               Degree = c(compare_strength$Before.HAB_degree, compare_strength$After.HAB_degree),
                               Strength = c(compare_strength$Before.HAB_strength, compare_strength$After.HAB_strength))

# Add HI_type column
local_metrics_HI_1 <- local_metrics_HI[local_metrics_HI$Period == "1-Before.HAB", ]
local_metrics_HI_2 <- local_metrics_HI[local_metrics_HI$Period == "2-After.HAB", ]

## Initialize a new dataframe to store the results
result_df_1 <- data.frame()
## Initialize a counter
counter <- 1
for (i in HI_list) {
  index <- local_metrics_HI_1$ID %in% i[[1]]
  # Create a new row for each ID that falls into four categories
  new_rows <- local_metrics_HI_1[index, ]
  new_rows$HI <- names(HI_list[counter])
  # Increment the counter
  counter <- counter + 1
  # Append the new rows to the result dataframe
  result_df_1 <- rbind(result_df_1, new_rows)
}

result_df_2 <- data.frame()
## Initialize a counter
counter <- 1
for (i in HI_list) {
  index <- local_metrics_HI_2$ID %in% i[[2]]
  # Create a new row for each ID that falls into four categories
  new_rows <- local_metrics_HI_2[index, ]
  new_rows$HI <- names(HI_list[counter])
  # Increment the counter
  counter <- counter + 1
  # Append the new rows to the result dataframe
  result_df_2 <- rbind(result_df_2, new_rows)
}

result_df <- rbind(result_df_1, result_df_2)

# Make period a binary variable and HI a factorial
result_df$Period <- as.numeric(as.factor(result_df$Period)) - 1
result_df$HI <- as.factor(result_df$HI)

# Save dataset
saveRDS(result_df, "result_df.RData")


###########################################################################
# PART 3: Run Model ---------------------------------------------

# Set up data
result_df <- readRDS("result_df.RData")

# Visualize data
result_df$Period <- ifelse(result_df$Period == 0, "Pre_HAB", "Post_HAB")
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
hist(c(scale(result_df$Between))) # continuous
hist(c(scale(result_df$Degree)))
hist(c(scale(result_df$Strength)))

# Make dummy variables
result_df$BG <- ifelse(result_df$HI == "BG", 1, 0)
result_df$FG <- ifelse(result_df$HI == "FG", 1, 0)
result_df$SD <- ifelse(result_df$HI == "SD", 1, 0)

# Look into nodal regression
## HI Behavior Combined Two Year Period ##
fit_mcmc.b <- MCMCglmm(Between ~ BG * Period + FG * Period + SD, data = result_df, nitt = 10000)
summary(fit_mcmc.b) # Might need to use a negative binomial dist

fit_mcmc.s <- MCMCglmm(Strength ~ BG * Period + FG * Period + SD, data = result_df, nitt = 10000)
summary(fit_mcmc.s)

fit_mcmc.d <- MCMCglmm(Degree ~ BG * Period + FG * Period + SD, data = result_df, nitt = 10000)
summary(fit_mcmc.d)

# Check for model convergence
model <- fit_mcmc.s
plot(model$Sol)
plot(model$VCV)

# Extract Posteriors
posterior <- model$Sol

# Plot the posterior distribution
mcmc_intervals(posterior, pars = c("(Intercept)", "BG", "FG", "SD",
                                   "Period", "BG:Period", "Period:FG"))
mcmc_areas(
  posterior, 
  pars = c("(Intercept)", "HIFG", "HINF", "HISD",
           "Period", "HIFG:Period", "HINF:Period", "HISD:Period"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)


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
