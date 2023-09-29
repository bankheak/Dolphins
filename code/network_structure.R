# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# SOCIAL NETWORK STRUCTURE
###########################################################################

# Set working directory here
setwd("../data")

# Add helpful functions
source("../code/functions.R") # edgelist function

###########################################################################
# PART 1: Structure Network ------------------------------------------------

## load all necessary packages
library(igraph) # Look at Dai Shizuka/Jordi Bascompte
library(tnet) # For weights
library(sna)
library(statnet)
library(doParallel)
library(png)

# Read in social association matrix
nxn <- readRDS("nxn.RData")
list_years <- readRDS("list_years.RData")
img.3 =readPNG("Dolphin.png") 

# Test one year at a time
year <- 1

## Create social network
ig <- lapply(nxn, function (df) {
  graph_from_adjacency_matrix(
  df,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE)})

# Set the node names based on row names
row_names <- lapply(nxn, function (df) {rownames(df)})
for (i in seq_along(ig)) {
  V(ig[[i]])$name <- row_names[[i]]
}

## Only show IDs of HI dolphins
### subset_HI in "GLMM.R"
HI_data <-  diff_raw(subset_HI(list_years))
row_names_HI <- lapply(HI_data, function (df) {
  as.vector(df$Code[(df$DiffHI == "BG" | df$DiffHI == "SD" | 
                                    df$DiffHI == "FG") & df$Freq > 0])})

# Plot network
# Set up the plotting area with 1 row and 2 columns for side-by-side plots
par(mfrow=c(1, 2))

# Loop through the list of graphs and plot them side by side
for (i in 1:length(ig)) {
  plot(ig[[i]],
       layout = layout_with_fr(ig[[i]]),
       edge.width = E(ig[[i]])$weight * 4,
       vertex.size = sqrt(igraph::strength(ig[[i]], vids = V(ig[[i]]), mode = c("all"), loops = TRUE) * 10),
       vertex.frame.color = NA,
       vertex.label.family = "Helvetica",
       vertex.label = ifelse(V(ig[[i]])$name %in% row_names_HI[[i]], V(ig[[i]])$name, NA),
       vertex.label.color = "black",
       vertex.label.cex = 0.8,
       vertex.label.dist = 2,
       vertex.frame.width = 0.01)
}

# rasterImage(img.3, xleft=0, xright=1.9, ybottom=0, ytop=1.5)

# Reset the plotting area to its default configuration
par(mfrow=c(1, 1))


###########################################################################
# PART 2: Network & Local Properties ------------------------------------------------

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
el <- readRDS("el_years.RData")

# Weighted clustering coefficients
cluster <- lapply(el, function (df) {clustering_local_w(df, measure=c("am", "gm", "mi", "ma", "bi"))})

## Betweenness centrality
between <- lapply(el, function (df) {betweenness_w(df, alpha=1)})

# Closeness centrality
close <- lapply(el, function (df) {closeness_w(df, alpha=1)})

# Degree and strength centrality
strength <- lapply(el, function (df) {degree_w(df, measure=c("degree","output"), type="out", alpha=1)})

###########################################################################
# PART 3: Network & Global Properties ------------------------------------------------

#' Breakdown: connectance = length(which(as.dist(orca_hwi)!=0))/(N*(N-1)/2)
#' Number of nodes (number of rows in the association matrix)
N = nrow(nxn[[year]])
#' Number of possible links: 
#' Nodes*(Nodes-1)/2: (-1 removes the node itself; /2 removes repetitions)
total = N*(N-1)/2
# Number of realized links: all non-zero cells in the association matrix
real = length(which(as.dist(nxn[[year]])!=0))
# Connectance: realized/total
real/total

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

## Modularity Q-value
modularity(dolphin_walk[[year]])
## Number of modules
groups(dolphin_walk[[year]])
## Membership of modules
membership(dolphin_walk[[year]])
## Save the edgelist into a new object
auxrand <- as.data.frame(el[[year]])

# Permutate the link weights
sample(auxrand$vw)
## Save in the auxrand object
auxrand[,3] <- sample(auxrand$vw)

# Calculate the modularity Q-value for a new permutated edge list
## Create a network from the list of nodes
igrand <- graph.edgelist(el[[year]][,1:2]) 
### Add link weights
E(igrand)$weight <- el[[year]][,2]
### Make undirected graph
igrand <- as.undirected(igrand)
## Permutate the link weights
E(igrand)$weight <- sample(E(igrand)$weight)
## Calculate modularity Q-value
rmod <- walktrap.community(igrand)
modularity(rmod)
## Number of modules
groups(rmod)
## Membership of modules
membership(rmod)

# Difference from our empirical data?
modularity(dolphin_walk[[year]])
modularity(rmod)

# Run modularity permutations 1000 times
iter = 1000
randmod = numeric()
for(i in 1:iter){
  # Save the edgelist into a new object
  auxrand <- el[[year]]
  # igraph format
  igrand <- graph.edgelist(auxrand[,1:2]) # Create a network from the list of nodes
  E(igrand)$weight <- auxrand[,2] # Add link weights
  igrand <- as.undirected(igrand) # Make undirected graph
  # Permutate the link weights
  E(igrand)$weight <- sample(E(igrand)$weight)
  # calculate the modularity Q-value
  rand_walk <- walktrap.community(igrand)
  randmod[i] <- modularity(rand_walk) # Save Q-value into a vector
}

## Calculate the 95% confidence interval (two-tailed test)
ci = quantile(randmod, probs=c(0.025, 0.975), type=2)

## Compare with the empirical Q-value
data.frame(Q=modularity(dolphin_walk[[year]]), LowCI=ci[1], HighCI=ci[2])

## Visualization random Q distribution
hist(randmod, xlim=c(0,0.6))
### Empirical Q-value
abline(v= modularity(dolphin_walk[[year]]), col="red")
### 2.5% CI
abline(v= ci[1], col="blue")
### 97.5% CI
abline(v= ci[2], col="blue")


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


# Set the node names based on row names
BG <- SD <- FG <- vector("list", length = length(dolp_ig))

for (i in seq_along(dolp_ig)) {
  # Set the node names
  V(dolp_ig[[i]])$name <- rownames(nxn[[i]])
  
  ## Parse out what HI behavior they engage in
  BG[[i]] <- as.vector(HI_data[[i]]$Code[HI_data[[i]]$DiffHI == "BG" & HI_data[[i]]$Freq > 0])
  SD[[i]] <- as.vector(HI_data[[i]]$Code[HI_data[[i]]$DiffHI == "SD" & HI_data[[i]]$Freq > 0])
  FG[[i]] <- as.vector(HI_data[[i]]$Code[HI_data[[i]]$DiffHI == "FG" & HI_data[[i]]$Freq > 0])
  
  ## Initialize label_color attribute for each node
  V(dolp_ig[[i]])$label_color <- "black"
  
  ## Make a different text color for each category
  V(dolp_ig[[i]])$label_color[V(dolp_ig[[i]])$name %in% BG[[i]]] <- "red"  
  V(dolp_ig[[i]])$label_color[V(dolp_ig[[i]])$name %in% SD[[i]]] <- "yellow" 
  V(dolp_ig[[i]])$label_color[V(dolp_ig[[i]])$name %in% FG[[i]]] <- "blue"
}
# BGSD <- intersect(BG, SD)
# BGFG <- intersect(BG, FG)
# SDFG <- intersect(SD, FG)
# BGSDFG <- intersect(BGSD, FG)
# V(dolp_ig[[year]])$label_color[V(dolp_ig[[year]])$name %in% BGSD] <- "orange"  
# V(dolp_ig[[year]])$label_color[V(dolp_ig[[year]])$name %in% BGFG] <- "purple" 
# V(dolp_ig[[year]])$label_color[V(dolp_ig[[year]])$name %in% SDFG] <- "green"  
# V(dolp_ig[[year]])$label_color[V(dolp_ig[[year]])$name %in% BGSDFG] <- "brown" 

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
par(mfrow=c(1, 2))
# Main labels for the plots
main_labels <- c("1993-2004 Network", "2005-2014 Network")  # Replace with appropriate main labels

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
