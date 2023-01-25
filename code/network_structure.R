# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# SOCIAL NETWORK STRUCTURE
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins")

###########################################################################
# PART 1: Applying NBDA to 'begging' data ---------------------------------

## load all necessary packages
require(igraph)

# Read in social association matrix
nxn<- read.csv("nxn.csv")

## Create social network
ig <- graph_from_adjacency_matrix(as.matrix(nxn),
                                  mode = c("undirected"),
                                  weighted = TRUE,
                                  diag = F,
                                  add.colnames = NULL,
                                  add.rownames = NA)

