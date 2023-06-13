# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# SOCIAL NETWORK STRUCTURE
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: Structure Network ------------------------------------------------

## load all necessary packages
require(igraph) # Look at Dai Shizuka/Jordi Bascompte
require(tnet) # For weights
require(sna)
require(statnet)
require(doParallel)

# Read in social association matrix
nxn <- readRDS("nxn.RData")
sample_data <- read.csv("sample_data.csv")

# Test one year at a time
year <- 1

## Create social network
ig <- graph_from_adjacency_matrix(as.matrix(nxn[[year]]),
                                  mode = c("undirected"),
                                  weighted = TRUE,
                                  diag = F, # No loops
                                  add.colnames = T,
                                  add.rownames = NA)

# Plot network
plot(ig,
     layout = layout_with_fr(ig),
     # link weight, rescaled for better visualization
     edge.width= E(ig)$weight*4,
     # node size as degree (rescaled)
     vertex.size= sqrt(igraph::strength(ig, vids = V(ig), mode = c("all"), loops = TRUE) *10 ),
     vertex.frame.color= NA, #"black",
     vertex.label.family = "Helvetica",
     vertex.label.color="black", 
     vertex.label.cex=0.8, 
     vertex.label.dist=2, 
     # edge.curved=0,
     vertex.frame.width=0.01,
)

###########################################################################
# PART 2: Network & Global Properties ------------------------------------------------

# Edgelist: Nodes (i & j) and edge (or link) weight
source("../code/functions.R") # SRI & null permutation
## Edgelist for each year
years <- unique(sample_data$Year)
n.cores <- detectCores()
system.time({
  registerDoParallel(n.cores)
  el_years <- list()
  for (i in 1:length(years)) {
    el_years[[i]] <- matrix_to_edgelist(nxn[[i]], rawdata = FALSE, idnodes = FALSE)
  }   
})

saveRDS(el_years, "el_years.RData")
el <- readRDS("../data/el_years.RData")

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

###########################################################################
# PART 3: Permutate Link Weights ------------------------------------------------

## igraph format with weight
system.time({
  registerDoParallel(n.cores)
  dolphin_ig <- list()
  for (j in 1:length(years)) {
    dolphin_ig[[j]] <- graph.adjacency(as.matrix(nxn[[j]]),mode="undirected",weighted=TRUE,diag=FALSE)
  }  
  ### End parallel processing
  stopImplicitCluster()
})

# Modularity by the WalkTrap algorithm 
system.time({
  registerDoParallel(n.cores)
  dolphin_walk <- list()
  for (k in 1:length(years)) {
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
## Link weight distribution
auxrand$vw

# Permutate the link weights
sample(auxrand$vw)
## Save in the auxrand object
auxrand[,3] <- sample(auxrand$vw)

# Calculate the modularity Q-value for a new permutated edge list
## Create a network from the list of nodes
igrand <- graph.edgelist(el[[year]][,1:2]) 
### Add link weights
E(igrand)$weight <- el[[year]][,3]
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
  E(igrand)$weight <- auxrand[,3] # Add link weights
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
hist(randmod, xlim=c(0.2,0.6))
### Empirical Q-value
abline(v= modularity(dolphin_walk[[year]]), col="red")
### 2.5% CI
abline(v= ci[1], col="blue")
### 97.5% CI
abline(v= ci[2], col="blue")


###########################################################################
# PART 4: Modularity ------------------------------------------------

# Create an unweighted network
system.time({
  registerDoParallel(n.cores)
  dolp_ig <- list()
  for (l in 1:length(years)) {
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
newman <- cluster_leading_eigen(dolp_ig[[year]], steps = -1, weights = E(dolp_ig[[year]])$weight, 
                                start = NULL, options = arpack_defaults, callback = NULL, 
                                extra = NULL, env = parent.frame())
 

# Random color scheme
col <- rgb(runif(max(newman$membership)), 
           runif(max(newman$membership)), 
           runif(max(newman$membership)))

# Assign a random color to individuals of each module ('module')
V(dolp_ig[[year]])$color <- NA
for (i in 1:max(newman$membership)){
  sample(col)
  V(dolp_ig[[year]])$color[which(newman$membership==i)] = col[i]
}

plot(dolp_ig[[year]])

