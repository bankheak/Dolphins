# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# MANTEL TESTS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: Structure Data for ILV ------------------------------------------------

## load all necessary packages
require(ade4) # Look at Dai Shizuka/Jordi Bascompte
require(ncf) # For weights

# Read in social association matrix
gbi<- read.csv("gbi.csv")
source("../code/functions.R") # SRI & null permutation
nxn<- SRI.func(gbi)
nxn<-as.matrix(nxn)

# Transforming SRI similarity into distance
dolp_dist = nxn + 0.00001
dolp_dist <- 1-nxn
## Remove the redundant cells and the diagonal 
dolp_dist <- as.dist(dolp_dist)

# Read file in to retain ILV
orig_data <- read.csv("firstgen_data.csv")
## Make date into a date class
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")
## Group each individual by date and sighting
area_vect <- cbind(orig_data[,c(8,17)]) # Seperate subarea and ID
area_vect <- rbind(area_vect[1:100,]) # To check that calculations are correct
area_vect <- transform(area_vect, Area = as.numeric(factor(Subarea)))
area_vect <- cbind(area_vect[,c(2,3)]) # Seperate subarea and Code
area_vect <- as.vector(area_vect$Area)
axi<- get_group_by_individual(area_vect, data_format = "individuals")

# Create simularity matrix
nxn_area<- sim.func(axi)

