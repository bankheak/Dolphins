# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# SOCIAL NETWORK STRUCTURE
###########################################################################

# Set working directory here
setwd("../data")

# Add helpful functions
source("../code/functions.R") # edgelist function and diff_raw(subset_HI())

###########################################################################
# PART 1: Structure Network ------------------------------------------------

## load all necessary packages
library(igraph) # Measure centrality here
library(tnet) # For weights
library(sna)
library(statnet)
library(doParallel) # For faster computing
library(ggplot2)
library(gridExtra) # To combine plots
library(reshape) # To rearrange a data frame
library(cowplot) # To add a legend

# Read in social association matrix
nxn <- readRDS("nxn.RData")
list_years <- readRDS("list_years.RData")

# Read in overlapped social association matrix
nxn_ovrlap <- readRDS("nxn_ovrlap.RData")
list_years_ovrlap <- readRDS("list_years_ovrlap.RData")

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
ig_ovrlap <- ig_func(nxn_ovrlap)

# Set the node names based on row names
row_name_assign <- function(nxn, ig) {
  row_names <- lapply(nxn, function (df) {rownames(df)})
  for (i in seq_along(ig)) {
    V(ig[[i]])$name <- row_names[[i]]
  }
}

row_name_assign(nxn, ig)
row_name_assign(nxn_ovrlap, ig_ovrlap)

## Only show IDs of HI dolphins
row_name_HI_func <- function(list_years) {
  HI_data <-  diff_raw(subset_HI(list_years))
  row_names_HI <- lapply(HI_data, function (df) {
    as.vector(df$Code[(df$DiffHI == "BG" | df$DiffHI == "SD" | 
                         df$DiffHI == "FG") & df$Freq > 0])})
  return(row_names_HI)
}

row_names_HI <- row_name_HI_func(list_years)
row_names_HI_ovrlap <- row_name_HI_func(list_years_ovrlap)

# Plot network
# Set up the plotting area with 1 row and 2 columns for side-by-side plots
par(mfrow=c(1, 2), mar = c(0.5, 0.5, 0.5, 0.5))

main_labels <- c("1993-2004 Network", "2005-2014 Network")

# Loop through the list of graphs and plot them side by side
for (i in 1:length(ig)) {
  plot(ig[[i]],
       layout = layout_with_fr(ig[[i]]),
       edge.width = E(ig[[i]])$weight * 4, # edge thickness
       vertex.size = sqrt(igraph::strength(ig[[i]], vids = V(ig[[i]]), mode = c("all"), loops = TRUE) * 10), # Changes node size based on an individuals strength (centrality)
       vertex.frame.color = NA,
       vertex.label.family = "Helvetica",
       vertex.label = ifelse(V(ig[[i]])$name %in% row_names_HI[[i]], V(ig[[i]])$name, NA),
       vertex.label.color = "black",
       vertex.label.cex = 0.8,
       vertex.label.dist = 2,
       vertex.frame.width = 0.01)
  # Add the main label above the plot
  title(main = main_labels[i], line = -1)
}

# rasterImage(img.3, xleft=0, xright=1.9, ybottom=0, ytop=1.5)

# Reset the plotting area to its default configuration
par(mfrow=c(1, 1))


###########################################################################
# PART 2: Network & Local Properties ------------------------------------------------

nxn <- nxn_ovrlap
list_years <- list_years_ovrlap
# Edgelist: Nodes (i & j) and edge (or link) weight
n.cores <- detectCores()
system.time({
  registerDoParallel(n.cores)
  el_years <- list()
  for (i in seq_along(list_years)) {
    el_years[[i]] <- matrix_to_edgelist(nxn[[i]], rawdata = FALSE, idnodes = FALSE)
  }  
  ### End parallel processing
  stopImplicitCluster()
})

saveRDS(el_years, "el_years.RData")
saveRDS(el_years, "el_years_ovrlap.RData")

el <- readRDS("el_years.RData")
el_ovrlap <- readRDS("el_years_ovrlap.RData")

el <- el_ovrlap
row_names_HI <- row_names_HI_ovrlap

# Set the node names based on row names
get_names <- function (matrix, metric) {
  row_names <- lapply(matrix, function (df) {rownames(df)})
for (i in seq_along(metric)) {
  metric[[i]][,1] <- row_names[[i]]
}
  return(metric)
  }

# Insert the transitioned IDs to the first period
transitioned_IDs <- setdiff(row_names_HI[[2]], row_names_HI[[1]])
all_HI_IDs <- c(transitioned_IDs, row_names_HI[[1]])

# Betweenness centrality
between <- lapply(el, function (df) {betweenness_w(df, alpha=1)})
between_diffs <- get_names(nxn, between)
between_diffs_HI <- lapply(seq_along(between_diffs), function(i) {
  df <- between_diffs[[i]]
  df_new <- as.data.frame(df[df[, 1] %in% all_HI_IDs, , drop = FALSE])
  return(df_new)
})
compare_between <- merge(
  between_diffs_HI[[1]], 
  between_diffs_HI[[2]], 
  by.x = "node", 
  by.y = "node"
)
colnames(compare_between) <- c("ID", "Period.1", "Period.2")
compare_between[, c(2, 3)] <- sapply(compare_between[, c(2, 3)], as.numeric)

# Degree and strength centrality
strength <- lapply(el, function (df) {degree_w(df, measure=c("degree","output"), type="out", alpha=1)})
strength_diffs <- get_names(nxn, strength)
strength_diffs_HI <- lapply(seq_along(strength_diffs), function(i) {
  df <- strength_diffs[[i]]
  df_new <- as.data.frame(df[df[, 1] %in% all_HI_IDs, , drop = FALSE])
  return(df_new)
})
compare_strength <- merge(
  strength_diffs_HI[[1]], 
  strength_diffs_HI[[2]], 
  by.x = "node", 
  by.y = "node"
)
colnames(compare_strength) <- c("ID", "Period.1_degree", "Period.1_strength", "Period.2_degree", "Period.2_strength")
compare_strength[, c(2:5)] <- sapply(compare_strength[, c(2:5)], as.numeric)

# Look at all of the local metrics together
HI_data <-  diff_raw(subset_HI(list_years))

## Add a column containing HI type
names_BG <- unlist(lapply(HI_data, function (df) {
  as.vector(df$Code[df$DiffHI == "BG" & df$Freq > 0])}))
names_SD <- unlist(lapply(HI_data, function (df) {
  as.vector(df$Code[df$DiffHI == "SD" & df$Freq > 0])}))
names_FG <- unlist(lapply(HI_data, function (df) {
  as.vector(df$Code[df$DiffHI == "FG" & df$Freq > 0])}))

# Combine the data
local_metrics_HI <- data.frame(ID = compare_between$ID,
                               Period = c(rep("Period.1", nrow(compare_between)), rep("Period.2", nrow(compare_between))),
                               Between = c(compare_between$Period.1, compare_between$Period.2),
                               Strength = c(compare_strength$Period.1_strength, compare_strength$Period.2_strength))

## Add a rown to compare the averages of each metric with HI IDs
avg_metrics <- data.frame(ID = "Average",
                          Period = c("Period.1", "Period.2"),
                          Between = c(mean(between[[2]][, 2]), mean(between[[1]][, 2])),
                          Strength = c(mean(strength[[2]][, 3]), mean(strength[[1]][, 3])))

local_metrics_HI <- rbind(local_metrics_HI, avg_metrics)

# Add HI_type column
local_metrics_HI$HI_type <- ifelse(local_metrics_HI$ID %in% names_BG, "BG", 
                            ifelse(local_metrics_HI$ID %in% names_SD, "SD", 
                            ifelse(local_metrics_HI$ID %in% names_FG, "FG", "NA")))

# Reshape the data from wide to long format
local_metrics_HI <- melt(local_metrics_HI, id.vars = c("ID", "HI_type", "Period"), variable.name = "Metric")
colnames(local_metrics_HI) <- c("ID", "HI_type", "Period", "Metric", "value")

# Make sure metric is in character
local_metrics_HI$Metric <- as.character(local_metrics_HI$Metric)

# Get rid of the average values
local_met_HI <- local_metrics_HI[local_metrics_HI$HI_type != "NA", ]

# Plot individual metrics
plot_list <- list()
unique_metrics <- unique(local_met_HI$Metric)

for (i in seq_along(unique_metrics)) {
  metric <- unique_metrics[i]
  
  # Filter data for the current metric
  metric_data <- local_met_HI[local_met_HI$Metric == metric,]
  
  # Get the corresponding value for NA, Period.1 and the current metric
  value_na_period1 <- local_metrics_HI$value[local_metrics_HI$HI_type == "NA" & 
                                               local_metrics_HI$Period == "Period.1" & 
                                               local_metrics_HI$Metric == metric]
  
  # Get the corresponding value for NA, Period.2 and the current metric
  value_na_period2 <- local_metrics_HI$value[local_metrics_HI$HI_type == "NA" & 
                                               local_metrics_HI$Period == "Period.2" & 
                                               local_metrics_HI$Metric == metric]
  
  # Create the plot with points
  current_plot <- ggplot(metric_data, aes(x = HI_type, y = value, color = Period, label = ID)) +
    geom_point(position = position_dodge(width = 0.5)) +  # Adjust width as needed
    labs(x = "HI Type", y = NULL, color = "Period") +
    ggtitle(paste(metric)) +
    theme(panel.background = element_blank()) +
    geom_hline(yintercept = value_na_period1, col = "red", linetype = "dashed") +
    geom_hline(yintercept = value_na_period2, col = "blue", linetype = "dashed") +
    theme(legend.position = "none") +
    geom_text(position = position_dodge(width = 0.5), vjust = -0.5)  # Adjust vjust as needed
  
  # Manually set colors for points
  current_plot <- current_plot + scale_color_manual(values = ifelse(metric_data$ID %in% transitioned_IDs, "black", ifelse(metric_data$Period == "Period.2", "blue", "red")))
  
  plot_list[[i]] <- current_plot
}

# Arrange plots side by side
grid.arrange(grobs = plot_list, ncol = 2)

###########################################################################
# PART 3: Network & Global Properties ------------------------------------------------

#' Breakdown: connectance = length(which(as.dist(orca_hwi)!=0))/(N*(N-1)/2)
## Calculate connectance for each matrix
calculate_connectance <- function(matrix) {
  N <- nrow(matrix)
  total <- N * (N - 1) / 2
  real <- sum(matrix != 0)  # Count non-zero elements
  connectance <- real / total
  return(connectance)
}

connectance_list <- lapply(nxn, calculate_connectance)

# Shortest path lengths (geodesics) and diameter
# # mean shortest path
dist <- lapply(ig, function(df) {mean_distance(df)})
dist

###########################################################################
# PART 4: Permutate Link Weights ------------------------------------------------

## igraph format with weight
system.time({
  registerDoParallel(n.cores)
  dolphin_ig <- list()
  for (j in seq_along(list_years)) {
    dolphin_ig[[j]] <- graph.adjacency(as.matrix(nxn[[j]]),
                                       mode="undirected",
                                       weighted=TRUE, diag=FALSE)
  }  
  ### End parallel processing
  stopImplicitCluster()
})

# Modularity by the WalkTrap algorithm 
system.time({
  registerDoParallel(n.cores)
  dolphin_walk <- list()
  for (k in seq_along(list_years)) {
    dolphin_walk[[k]] <- cluster_walktrap(dolphin_ig[[k]], weights = E(dolphin_ig[[k]])$weight, 
                                     steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)
  } 
  ### End parallel processing
  stopImplicitCluster()
})

# Run modularity permutations 1000 times for each matrix
run_mod <- function(el, dolphin_walk_list) {
  iter <- 1000
  randmod <- numeric(iter)  # Initialize a numeric vector to store Q-values
  
  for (i in 1:iter) {
    # Save the edgelist into a new object and permutate the link weights
    auxrand <- el
    auxrand[, 2] <- sample(auxrand[, 2])
    
    # Create an igraph graph from the permuted edgelist
    igrand <- graph.edgelist(as.matrix(auxrand[, 1:2]), directed = FALSE)
    E(igrand)$weight <- auxrand[, 2]  # Assign link weights
    
    # Calculate modularity using walktrap community detection
    rand_walk <- walktrap.community(igrand)
    randmod[i] <- modularity(rand_walk)  # Save Q-value into the vector
  }
  
  # Calculate the 95% confidence interval (two-tailed test)
  ci <- quantile(randmod, probs = c(0.025, 0.975), type = 2)
  
  # Visualization of the random Q distribution
  hist(randmod, xlim = c(0, 0.6), main = "Random Q Distribution", xlab = "Q-value", ylab = "Frequency", col = "lightblue")
  
  # Empirical Q-value
  abline(v = modularity(dolphin_walk_list), col = "red")
  
  # 2.5% CI
  abline(v = ci[1], col = "blue")
  
  # 97.5% CI
  abline(v = ci[2], col = "blue")
  
  # Return a data frame with Q-value and confidence intervals
  result <- data.frame(Q = modularity(dolphin_walk_list), LowCI = ci[1], HighCI = ci[2])
  return(result)
}

run_mod(el_list = el[[1]], dolphin_walk_list = dolphin_walk[[1]])
run_mod(el_list = el[[2]], dolphin_walk_list = dolphin_walk[[2]])


###########################################################################
# PART 5: Modularity ------------------------------------------------

# Create an unweighted network
system.time({
  registerDoParallel(n.cores)
  dolp_ig <- list()
  for (l in seq_along(list_years)) {
    dolp_ig[[l]] <- graph.edgelist(el[[l]][,1:2])
    # Add the edge weights to this network
    E(dolp_ig[[l]])$weight <- as.numeric(el[[l]][,3])
    # Create undirected network
    dolp_ig[[l]] <- as.undirected(dolp_ig[[l]])
  }   
  ### End parallel processing
  stopImplicitCluster()
})

# Newman's Q modularity
newman <- lapply(dolp_ig, function (df) {cluster_leading_eigen(df, steps = -1, weights = E(df)$weight, 
                                start = NULL, options = arpack_defaults, callback = NULL, 
                                extra = NULL, env = parent.frame())})

# Set the node names and label colors based on HI behavior
BG <- SD <- FG <- BGSD <- BGFG <- SDFG <- BGSDFG <- vector("list", length = length(dolp_ig))

for (i in seq_along(dolp_ig)) {
  # Set the node names
  V(dolp_ig[[i]])$name <- rownames(nxn[[i]])
  
  # Parse out what HI behavior they engage in
  BG[[i]] <- as.vector(HI_data[[i]]$Code[HI_data[[i]]$DiffHI == "BG" & HI_data[[i]]$Freq > 0])
  SD[[i]] <- as.vector(HI_data[[i]]$Code[HI_data[[i]]$DiffHI == "SD" & HI_data[[i]]$Freq > 0])
  FG[[i]] <- as.vector(HI_data[[i]]$Code[HI_data[[i]]$DiffHI == "FG" & HI_data[[i]]$Freq > 0])
  BGSD[[i]] <- intersect(BG[[i]], SD[[i]])
  BGFG[[i]] <- intersect(BG[[i]], FG[[i]])
  SDFG[[i]] <- intersect(SD[[i]], FG[[i]])
  BGSDFG[[i]] <- intersect(BGSD[[i]], FG[[i]])
  
  # Initialize label_color attribute for each node
  V(dolp_ig[[i]])$label_color <- "black"
  
  # Set label colors based on categories
  node_names <- V(dolp_ig[[i]])$name
  V(dolp_ig[[i]])$label_color <- ifelse(node_names %in% BGSDFG[[i]], "brown",
                                        ifelse(node_names %in% BGFG[[i]], "purple",
                                               ifelse(node_names %in% SDFG[[i]], "green",
                                                      ifelse(node_names %in% BGSD[[i]], "orange",
                                                             ifelse(node_names %in% FG[[i]], "blue",
                                                                    ifelse(node_names %in% SD[[i]], "yellow",
                                                                           ifelse(node_names %in% BG[[i]], "red", "black")))))))
  
}

# Generate a vector of colors based on the number of unique memberships
for (i in seq_along(dolp_ig)) {
  V(dolp_ig[[i]])$color <- NA
  col <- rainbow(max(newman[[i]]$membership))
  
  for (j in 1:max(newman[[i]]$membership)){
    V(dolp_ig[[i]])$color[which(newman[[i]]$membership==j)] <- col[j]
  }
}

# Make sure the HI dolphins stand out
for (i in seq_along(dolp_ig)) {
  V(dolp_ig[[i]])$size <- ifelse(V(dolp_ig[[i]])$name %in% row_names_HI[[i]], 10, 5)
}

# Set up the plotting area with 1 row and 2 columns for side-by-side plots
par(mfrow=c(1, 2), mar = c(0.5, 0.5, 2, 0.5))

# Main labels for the plots
main_labels <- c("1993-2004 Network", "2005-2014 Network")

# Plot the graph with individual IDs as labels
for (i in seq_along(dolp_ig)) {
plot(dolp_ig[[i]],
     layout = layout_with_fr(dolp_ig[[i]]),
     # link weight, rescaled for better visualization
     edge.width= E(dolp_ig[[i]])$weight*4,
     # node size as degree (rescaled)
     vertex.size= V(dolp_ig[[i]])$size,
     vertex.frame.color= NA, #"black",
     vertex.label.family = "Helvetica",
     vertex.label=ifelse(V(dolp_ig[[i]])$name %in% row_names_HI[[i]], as.character(V(dolp_ig[[i]])$name), NA), 
     vertex.label.color = V(dolp_ig[[i]])$label_color, 
     vertex.label.cex=0.8, 
     vertex.label.dist=0.5, 
     # edge.curved=0,
     vertex.frame.width=0.01)
  # Add the main label above the plot
  title(main = main_labels[i], line = -1)
  }
