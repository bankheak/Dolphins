# Network Analysis of Anthropogenic Influences in Bottlenose dolphins

# Local Network Analysis Hypothesis #2 #

# Set working directory here
setwd("../data")

# Load all necessary packages
library(tnet) # For weights
library(igraph) # Measure centrality here
library(ggraph)
library(grid)
library(assortnet) # associative indices
library(ggplot2) # Visualization
library(abind) # array
library(MCMCglmm) # MCMC models
library(coda)
library(bayesplot) # plot parameters
library(doParallel)
library(hrbrthemes) # plot themes
library(viridis) # plot themes
library(patchwork) # plotting together
library(car) # durbinWatsonTest
library(DescTools) #Schaff post-hoc test
library(lme4) # lmm
library(nlme) # unequal variance weights
library(lmerTest) # summary output of lmm
library(emmeans) # post-hoc test
library(effects) # visualize effects
library(sjPlot) # Confidence intervals
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

## Create social network
ig_func <- function(nxn) {
  ig <- lapply(nxn, function (df) {
    graph_from_adjacency_matrix(
      df,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE)})
  return(ig)}

ig <- ig_func(nxn)

# Set the node names based on row names
row_name_assign <- function(nxn, ig) {
  row_names <- lapply(nxn, function (df) {rownames(df)})
  for (i in seq_along(ig)) {
    V(ig[[i]])$name <- row_names[[i]]
  }
}

row_name_assign(nxn, ig)

# Save ig object
saveRDS(ig, "ig.RData")

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

#Eigen centrality
eigen <- lapply(ig, function (df) {eigen_centrality(df)})

eigen <- lapply(eigen, function (df) {
  eigen <- as.data.frame(df$vector)
  eigen$ID <- rownames(eigen)
  colnames(eigen) <- c("eigen_cent", "ID")
  return(eigen)})

compare_eigen <- merge(
  merge(eigen[[1]], eigen[[2]], by = "ID"),
  eigen[[3]], by = "ID"
)

colnames(compare_eigen) <- c("ID", "Before_HAB", "During_HAB", "After_HAB")
compare_eigen[, c(2:4)] <- sapply(compare_eigen[, c(2:4)], as.numeric)

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

saveRDS(HI_list, "HI_list.RData")

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
                                            compare_strength$After_HAB_strength),
                               Eigen = c(compare_eigen$Before_HAB, 
                                         compare_eigen$During_HAB,
                                         compare_eigen$After_HAB))

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

# Create weighted centrality metric

## Scale the centrality metrics
result_df$Degree <- c(scale(result_df$Degree))
result_df$Strength <- c(scale(result_df$Strength))
result_df$Between <- c(scale(result_df$Between))
result_df$Eigen <- c(scale(result_df$Eigen))

## Put weights on centrality
weighted_degree <- result_df$Degree *  0.25
weighted_strength <- result_df$Strength *  0.25
weighted_betweenness <- result_df$Between *  0.25
weighted_eigen <- result_df$Eigen * 0.25

##Create the weighted metric
result_df$composite_centrality <- weighted_degree + weighted_strength + 
                                  weighted_betweenness + weighted_eigen

# Save dataset
saveRDS(result_df, "result_df.RData")


###########################################################################
# PART 3: Run Model ---------------------------------------------

# Read in data
result_df <- readRDS("result_df.RData")

# First Visualize data
ggplot(result_df, aes(x = Period, y = composite_centrality, fill = HI)) + 
  geom_boxplot() 

# Make dummy variables
result_df$BG <- ifelse(result_df$HI == "BG", 1, 0)
result_df$FG <- ifelse(result_df$HI == "FG", 1, 0)
result_df$SD <- ifelse(result_df$HI == "SD", 1, 0)
# Make factor variables
result_df$HI <- as.factor(result_df$HI)
result_df$Period <- as.factor(result_df$Period)

# Check assumptions of model
test_model <- lm(composite_centrality ~ BG * During + BG * After +
                   FG * During + FG * After + 
                   SD * During + SD * After, data = result_df)
summary(test_model)
## Check distributions
hist(result_df$composite_centrality) # normal
## Check for variance among groups
bartlett.test(composite_centrality ~ HI, data = result_df) # not equal
## Independent
durbinWatsonTest(test_model) # not independent

# Make ID numeric
result_df$numeric_ID <- as.numeric(factor(result_df$ID))

# Fit the LMM
lmm_model_0 <- lme(composite_centrality ~  1, random = ~1 | numeric_ID,
                   weights = varIdent(form = ~1 | HI), data = result_df)
lmm_model_1 <- lme(composite_centrality ~  BG + FG + SD, 
                   random = ~1 | numeric_ID, weights = varIdent(form = ~1 | HI), 
                   data = result_df)
lmm_model_2 <- lme(composite_centrality ~  BG + FG + SD + During + After, 
                   random = ~1 | numeric_ID, weights = varIdent(form = ~1 | HI), 
                   data = result_df)
lmm_model_3 <- lme(composite_centrality ~  BG * During + BG * After + 
                     FG * During + FG * After + SD * During + SD *After, 
                   random = ~1 | numeric_ID, weights = varIdent(form = ~1 | HI), 
                   data = result_df)

# Fit the MCMC
fit_mcmc_0 <- MCMCglmm(composite_centrality ~ 1, random=~ numeric_ID, rcov = ~ HI,
                       data = result_df, nitt = 20000) 
fit_mcmc_1 <- MCMCglmm(composite_centrality ~ BG + FG + SD, 
                       random=~ numeric_ID, rcov = ~ HI, data = result_df, nitt = 20000) 
fit_mcmc_2 <- MCMCglmm(composite_centrality ~ BG + FG + SD + During + After, 
                        random=~ numeric_ID, rcov = ~ HI, data = result_df, nitt = 20000) 
fit_mcmc_3 <- MCMCglmm(composite_centrality ~ BG * During + BG * After + 
                         FG * During + FG * After + SD * During + SD * After, 
                       random=~ numeric_ID, rcov = ~ HI, data = result_df, nitt = 20000) 

summary(fit_mcmc_0)
summary(fit_mcmc_1) 
summary(fit_mcmc_2) # Lowest DIC
summary(fit_mcmc_3)

# Model Selection
AIC(lmm_model_0, lmm_model_1, lmm_model_2, lmm_model_3)

# Print the summary of the model
summary(lmm_model_2)

# Run post-hoc test
emm_pairs <- emmeans(lmm_model_2, pairwise ~ BG + FG + SD + During + After, adjust = "hochberg")
summary(emm_pairs, infer = TRUE)

# Visualize effects
effects_lmm_model <- allEffects(lmm_model_2)
plot(effects_lmm_model)
plot_model(lmm_model_2)


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

###########################################################################
# PART 5: Multinetwork ---------------------------------------------

# Only show IDs of HI dolphins
HI_list <- readRDS("HI_list.RData")
HI_list <- HI_list[-4] # Get rid of natural foragers

# Read in ig object
ig <- readRDS("ig.RData")

# ---Plot network---
# Set up the plotting area with 1 row and 2 columns for side-by-side plots
# Initialize a list to store layout information for each graph
layout_list <- list()

# Loop through the list of graphs and save layout information
for (i in 1:length(ig)) {
  layout_list[[i]] <- layout_with_fr(ig[[i]])
}

# Set up the plotting layout
layout.matrix <- matrix(c(1:9), nrow = 3, ncol = 3)
layout(mat = layout.matrix)    
par(mar = c(0.6, 0.6, 0.6, 0.6))

# Loop through the list of graphs and plot them side by side
for (j in 1:length(HI_list)) {  # Loop through columns first
  
  # Extract layout for this graph
  combined_layout <- layout_list[[1]]
  counter <- 0
  
  for (i in 1:length(ig)) {  # Loop through rows
    
    counter <- counter + 1
    
    # Get nodes for each behavior
    labeled_nodes <- V(ig[[i]])$name %in% HI_list[[j]][[i]]  # Fixed index here
    
    # Create the plot
    plot(ig[[i]],
         layout = combined_layout,
         edge.width = E(ig[[i]])$weight * 4, # edge thickness
         edge.color = adjustcolor("grey", alpha.f = 0.2),
         vertex.size = ifelse(labeled_nodes, 10, 3), #sqrt(igraph::strength(ig[[i]], vids = V(ig[[i]]), mode = c("all"), loops = TRUE) * 10), # Changes node size based on an individual's strength (centrality)
         vertex.frame.color = NA,
         vertex.label.family = "Helvetica",
         vertex.label = ifelse(labeled_nodes, V(ig[[i]])$name, NA),
         vertex.label.color = "black",
         vertex.label.cex = 0.8,
         vertex.label.dist = 2,
         vertex.frame.width = 0.01,
         vertex.color = "black")
    
    # Add the plot with a box around it
    box()
    
  }
}


# Set up data
result_df <- readRDS("result_df.RData")

# Plot the density plots for each period
plots_list <- list()

for (i in 1:length(unique(result_df$Period))) {
  
  period_to_plot <- unique(result_df$Period)[i] # each period
  filtered_df <- subset(result_df, Period == period_to_plot) # Separate data
  
  mean_nf <- mean(filtered_df$composite_centrality[filtered_df$HI == "NF"], na.rm = TRUE) # Calculate mean for HI=="NF"
  
  plot <- ggplot(filtered_df[filtered_df$HI != "NF", ], aes(x = HI, y = composite_centrality, fill = HI)) +
    geom_violin(trim = FALSE, alpha = 0.4) + # Create violin plot
    geom_boxplot(width=0.1, color="black", alpha=0.2) +
    geom_jitter(width = 0.1, alpha = 0.6) + # Add jittered points for visibility
    geom_hline(yintercept = mean_nf, linetype = "dashed", color = "black", linewidth = 1.5) + # Add horizontal line
    theme_ipsum() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank()) 
  
  plots_list[[i]] <- plot
}


# Output plots
plots_list[[1]]
plots_list[[2]]
plots_list[[3]]


# Plot the density plots for each HI
# Define shades of red
red_shades <- c("#FFCCCC", "#FF0000", "#FF6666")

plots_list_HI <- list()

for (i in 1:(length(unique(result_df$HI))-1)) {
  
  HI_to_plot <- unique(result_df$HI)[i] # each HI category
  filtered_df <- subset(result_df, HI == HI_to_plot) # Separate data
  
  plot <- ggplot(filtered_df, aes(x = Period, y = composite_centrality, fill = as.factor(Period))) +
    geom_violin(trim = FALSE, alpha = 0.4) + # Create violin plot
    geom_boxplot(width = 0.1, color = "black", alpha = 0.2) +
    geom_jitter(width = 0.1, alpha = 0.6) + # Add jittered points for visibility
    scale_fill_manual(values = red_shades) + # Use manual scale for shades of red
    theme_ipsum() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    guides(fill = guide_legend(title = "Period")) # Add legend for Period
  
  plots_list_HI[[i]] <- plot
}

# Output plots
plots_list_HI[[1]]
plots_list_HI[[2]]
plots_list_HI[[3]]
