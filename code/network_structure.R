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

# Read in social association matrix
nxn<- read.csv("nxn.csv")
nxn<- as.matrix(nxn)

## Create social network
ig <- graph_from_adjacency_matrix(as.matrix(nxn),
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

# Source edgelist function
source("../code/functions.R")

# Edgelist: Nodes (i & j) and edge (or link) weight
el <- matrix_to_edgelist(nxn, rawdata = FALSE, idnodes = FALSE)

# Centrality measures
# Weighted clustering coefficients
clustering_local_w(el, measure=c("am", "gm", "mi", "ma", "bi"))

## Betweenness centrality
betweenness_w(el, alpha=1)       

# Closeness centrality
closeness_w(el, alpha=1)

# Degree and strength centrality
degree_w(el, measure=c("degree","output"), type="out", alpha=1)

# Connectance or Density: proportion of realized links
edge_density(ig, loops = FALSE)

#' Breakdown: connectance = length(which(as.dist(orca_hwi)!=0))/(N*(N-1)/2)
#' Number of nodes (number of rows in the association matrix)
N = nrow(nxn)
#' Number of possible links: 
#' Nodes*(Nodes-1)/2: (-1 removes the node itself; /2 removes repetitions)
total = N*(N-1)/2
# Number of realized links: all non-zero cells in the association matrix
real = length(which(as.dist(nxn)!=0))
# Connectance: realized/total
real/total

# Shortest path lengths (geodesics) and diameter: binary and weighted
# all binary shortest path lengths between nodes
distances(ig)

# mean shortest path
mean_distance(ig)

# Binary shortest path lengths
distance_table(ig, directed=FALSE)

# All Weighted Shortest path lenghts (or geodesics)
distance_w(el)
# weighted diameter
max(as.dist(distance_w(el)))

# Clustering Coefficient
# Weighted clustering coefficients (using 5 different methods). 
clustering_w (el, measure=c("am", "gm", "mi", "ma", "bi"))


###########################################################################
# PART 3: Modularity ------------------------------------------------

# Create an unweighted network
dolp_ig <- graph.edgelist(el[,1:2])
# Add the edge weights to this network
E(dolp_ig)$weight <- as.numeric(el[,3])
# Create undirected network
dolp_ig <- as.undirected(dolp_ig)

# Plot
plot(dolp_ig, edge.width=E(dolp_ig)$weight*4, vertex.size=10, vertex.label=NA, edge.curved=F)

# Newman's Q modularity
newman <- cluster_leading_eigen(dolp_ig, steps = -1, weights = E(dolp_ig)$weight, 
                                start = NULL, options = arpack_defaults, callback = NULL, 
                                extra = NULL, env = parent.frame())
plot(dolp_ig)

# Random color scheme
col <- rgb(runif(10), runif(10), runif(10))

# Assign a random color to individuals of each module ('module')
newman$membership

for (i in 1:length(newman$membership)){
  sample(col)
  V(dolp_ig)$color[which(newman$membership==i)] = col[i]
}
plot(dolp_ig)
