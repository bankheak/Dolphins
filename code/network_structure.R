# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# SOCIAL NETWORK STRUCTURE
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: Applying NBDA to 'begging' data ---------------------------------

## load all necessary packages
require(igraph) # Look at Dai Shizuka/Jordi Bascompte

# Read in social association matrix
nxn<- read.csv("nxn.csv")

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
