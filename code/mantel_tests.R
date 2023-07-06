# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# MANTEL TESTS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: Create HI Similarity Matrix  ------------------------------------------------

## load all necessary packages
require(ade4) # Look at Dai Shizuka/Jordi Bascompte
require(ncf) # For weights
require(vegan)

# Read file in to retain ILV
sample_data <- read.csv("sample_data.csv")

# Read in social association matrix and data
nxn <- readRDS("nxn.RData")
list_years <- readRDS("list_years.RData")

# Transforming SRI similarity into distance
year <- 1
dolp_dist = nxn[[year]] + 0.00001
dolp_dist <- 1-nxn[[year]]
## Remove the redundant cells and the diagonal 
dolp_dist <- as.dist(dolp_dist)

# Select variables from the raw data
aux <- sample_data[1:nrow(list_years[[year]]), 
                c('Code', 'Behaviors', 'HumanInteraction', 'ConfHI')]

# Use 'Behaviors' variable to extract "Feed" and create another variable with two classes (Feed, Other)
aux$Foraging <- "Other"
aux$Foraging[grepl(pattern = 'Feed',
                   x = aux$Behaviors,
                   ignore.case = FALSE, perl = FALSE,
                   fixed = FALSE, useBytes = FALSE)] = "Feed"
aux <- subset(aux, aux$Foraging == "Feed")
aux$ConfHI <- ifelse(aux$ConfHI == "0", 0, 1)

# Categorize ID to Foraging
IDbehav <- table(aux$Code, aux$Foraging)
IDbehav <- data.frame(IDbehav)[,c(1,3)]
colnames(IDbehav) <- c("Code", "Forg_Freq")

# Categorize ConfHI to IDs
rawHI <- as.matrix(table(aux$Code, aux$ConfHI))
rawHI <- data.frame(rawHI)
colnames(rawHI) <- c("Code", "ConfHI", "Freq")

## Add up the # of times each ID was seen in HI
IDbehav$HI <- rawHI$Freq[rawHI$ConfHI != 0]
IDdata <- IDbehav
colnames(IDdata) <- c("Code", "Foraging", "HI")

## Proportion of time FOraging spent in HI
IDdata$HIprop <- as.numeric(IDdata$HI)/as.numeric(IDdata$Foraging)

# Only ID to prop
HIprop_ID <- na.omit(IDdata[,c(1, 4)])

# ----Or I could measure HI per obs for each ID-----
IDdata_2 <- as.data.frame(table(aux$Code)) # FB07 (Group 4) & FB92 (Group 4) in year 1
IDdata_2$HI <- rawHI$Freq[rawHI[,2] != "0"]
IDdata_2$HIprop <- IDdata_2$HI/IDdata_2$Freq

# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
fake_HIprop <- HIprop_ID$HIprop
dissimilarity <- as.matrix(vegdist(fake_HIprop, method = 'euclidean'))
