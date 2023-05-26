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
require(vegan)

# Read file in to retain ILV
orig_data <- read.csv("sample_data.csv")

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
aux <- orig_data[1:nrow(list_years[[year]]), 
                c('Code', 'Behaviors', 'HumanInteraction', 'ConfHI')]

# Use 'Behaviors' variable to extract "Feed" and create another variable with two classes (Feed, Other)
aux$Foraging <- "Other"
aux$Foraging[grepl(pattern = 'Feed',
                   x = aux$Behaviors,
                   ignore.case = FALSE, perl = FALSE,
                   fixed = FALSE, useBytes = FALSE)] = "Feed"

# Categorize ID to Foraging
IDbehav <- table(aux$Code, aux$Foraging)
IDbehav <- data.frame(IDbehav)

# Categorize ConfHI to IDs
rawHI <- as.matrix(table(aux$Code, aux$ConfHI)[,1:2])
rawHI <- data.frame(rawHI)

# Take out the number of foraging events per ID
IDdata <- data.frame(Foraging = subset(IDbehav, IDbehav$Var2 == "Feed"))
IDdata <- IDdata[,c(1,3:4)]
## Add up the # of times each ID was seen in HI
IDdata$HI = rawHI$Freq[rawHI$Var2 == "P"]
## Proportion of time FOraging spent in HI
IDdata$HIprop = IDdata[,3]/IDdata[,2]

# How many observations do each ID have?
IDdata = as.data.frame(table(aux$Code))
IDdata$HI = as.vector(rowSums(rawHI))
IDdata$HIprop = IDdata$HI/IDdata$Freq

identical(row.names(rawHI), IDdata$Var1)


# Count how many times each individual was 'Feed'
# Count how many times each ind was in HI or not
# Divide the two to create the 'propHI'


aux$Behaviors[grepl(pattern = 'Feed',
                    x = aux$Behaviors,
                    ignore.case = FALSE, perl = FALSE,
                    fixed = FALSE, useBytes = FALSE)]

as.vector(table(aux$Behaviors))
table(aux$HumanInteraction)


# pick the feeding events
aux[which(aux$Behaviors == 'pFeed'), ]

# 
aux[which(aux$HumanInteraction != 999), ]




# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
fake_HIprop = t(varespec)[,1]
as.matrix(vegdist(fake_HIprop, method = 'euclidean'))


# Create simularity matrix
nxn_area<- sim.func(aux)

