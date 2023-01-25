# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# DYATIC SOCIAL ASSOCIATIONS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: Gambit of the group ---------------------------------------------

# Load all necessary packages
require(vegan)
require(asnipe)
require(assocInd)

# Read file in
group_data <- read.csv("93-14_data.csv")

# Make date into a date class
group_data$Date <- as.Date(as.character(group_data$Date), format="%d-%b-%y")

# Group each individual by date and sighting
group_data <- cbind(group_data[,c(2,11,17)])
group_data$Group <- cumsum(!duplicated(group_data[1:2]))
group_data <- cbind(group_data[,3:4])
sample_data <- cbind(group_data[,3:4]) # To check that calculations are correct
sample_data<- rbind(sample_data[1:10,])

# Gambit of the group
gbi<- get_group_by_individual(sample_data, data_format = "individuals")

# Create association index
source("../code/functions.R")
nxn<- SRI.func(gbi)
write.csv(nxn, "nxn.csv")
