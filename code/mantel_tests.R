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
kov <- readRDS("kov.RDS")

# Read in social association matrix and data
nxn <- readRDS("nxn.RData")
list_years <- readRDS("list_years.RData")

# Transforming SRI similarity into distance
year <- 5
dolp_dist = nxn[[year]] + 0.00001
dolp_dist <- 1-nxn[[year]]
## Remove the redundant cells and the diagonal 
dolp_dist <- as.dist(dolp_dist)

# Select variables from the raw data
data <- list_years[[year]]
aux <- data[, c('Code', 'Behaviors', 'HumanInteraction', 'ConfHI')]

length(unique(aux$Code)) # 381 individuals should stay consistent

# Use 'Behaviors' variable to extract "Feed" and create another variable with two classes (Feed, Other)
aux$Foraging <- "Other"
aux$Foraging[grepl(pattern = 'Feed',
                   x = aux$Behaviors,
                   ignore.case = FALSE, perl = FALSE,
                   fixed = FALSE, useBytes = FALSE)] = "Feed"
#aux <- subset(aux, aux$Foraging == "Feed")
aux$ConfHI <- ifelse(aux$ConfHI == "0", 0, 1)

# Categorize ID to Foraging
IDbehav <- table(aux$Code, aux$Foraging)
IDbehav <- as.data.frame(IDbehav, stringsAsFactors = FALSE)
IDbehav <- IDbehav[,c(1,3)]
colnames(IDbehav) <- c("Code", "Forg_Freq")
# Group by the 'Code' column and sum the frequencies
IDbehav <- aggregate(. ~ Code, data = IDbehav, sum)

# Categorize ConfHI to IDs
rawHI <- as.matrix(table(aux$Code, aux$ConfHI))
rawHI <- as.data.frame(rawHI, stringsAsFactors = FALSE)
colnames(rawHI) <- c("Code", "ConfHI", "Freq")

#HI <- subset(rawHI, subset=c(rawHI$ConfHI != 0))
#HI <- subset(HI, subset=c(HI$Freq != 0))

## Add up the # of times each ID was seen in HI
IDbehav$HI <- rawHI$Freq[rawHI$ConfHI != 0]
IDdata <- IDbehav
colnames(IDdata) <- c("Code", "Foraging", "HI")

## Proportion of time Foraging spent in HI
IDdata$HIprop <- as.numeric(IDdata$HI)/as.numeric(IDdata$Foraging)
IDdata[is.na(IDdata)] <- 0

# Only ID to prop
HIprop_ID <- IDdata[,c(1, 4)]

# ----Or I could measure HI per obs for each ID-----
#IDdata_2 <- as.data.frame(table(aux$Code)) # FB07 (Group 4) & FB92 (Group 4) in year 1
#IDdata_2$HI <- rawHI$Freq[rawHI[,2] != "0"]
#IDdata_2$HIprop <- IDdata_2$HI/IDdata_2$Freq

# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
fake_HIprop <- HIprop_ID$HIprop
dissimilarity_HI <- as.matrix(dist(as.matrix(fake_HIprop), method = "euclidean"))
dissimilarity_HI[is.na(dissimilarity_HI)] <- 0

dissimilarity_HI <- as.dist(dissimilarity_HI) # HI dissimilarity
kov <- as.dist(kov) # Home range overlap

# Make sure all matrices have the same dimensions
length(unique(list_years[[year]]$Code)) == length(unique(HIprop_ID$Code))


# Dissimilarity matrix
mantel.rtest(kov, dolp_dist, nrepet=999)

mantel.rtest(dissimilarity_HI, dolp_dist, nrepet=999) 
