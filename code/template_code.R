# 'Multi-network Network-Based Diffusion Analysis

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins")

###########################################################################
# PART 1: calculate dyadic home range overlaps ----------------------------

## load all necessary packages
require(rmarkdown)  # “Knit” button (Ctrl+Shift+K) displays preview
require(ggplot2)  # Plot designs


###########################################################################
# PART 2: applying NBDA to 'begging' data ---------------------------------

## load all necessary packages
install.packages("devtools")
require(devtools)
install_github("whoppitt/NBDA")
require(NBDA)

## Read in all network data and turn into matrix
### Vertical network
vert_all <- read.csv("social_vertical.csv", row.names=1, header=TRUE)
vert_all <- as.matrix(vert_all)

### Horizontal network
hor_all <- read.csv("social_horizontal.csv", row.names=1, header=TRUE)
hor_all <- as.matrix(hor_all)

### Ecological network
ecol_all <- read.csv("HR_overlaps.csv", row.names=1, header=TRUE)
ecol_all <- as.matrix(ecol_all)

### Relatedness network
relate_all <- read.csv("relatedness.csv", row.names=1, header=TRUE)
relate_all <- as.matrix(relate_all)

### ILVs
ILV_all <- read.csv("ILVs.csv", header=TRUE, sep=",")
ILV_all[ILV_all=="<NA>"]=NA

## Make list of IDs that have been seen at least 7 times
ILV <- subset(ILV_all, subset=ILV_all$Number_sightings>6) # at least 7 sightings
ILV <- subset(ILV, subset=ILV$Not_weaned==0) # separate weaned

IDs <- ILV$id_individual
length(IDs) #___ individuals remain

## Reduce vertical social network to only include the 415 individuals 
## that made the cut-off point
num <- which(colnames(vert_all) %in% IDs)
vert_all <- vert_all[num, num]
dim(vert_all)
class(vert_all)

  ## Repeat for horizontal social network
num <- which(colnames(hor_all) %in% IDs)
hor_all <- hor_all[num, num]
dim(hor_all)
class(hor_all)

  ## Repeat for ecological network
num <- which(colnames(ecol_all) %in% IDs)
ecol <- ecol_all[num, num]
dim(ecol)
class(ecol)

  ## Repeat for genetic network
num <- which(colnames(relate_all) %in% IDs)
relate <- relate_all[num, num]
dim(relate)
class(relate)

## Extract beggars (learners and demonstrators)
beggars <- subset(ILV, subset=ILV$Beggar=="yes")
beggars <- spongers[order(beggars$Order_acquisition),]

### ID codes of all beggars
beggar_all <- spongers$id_individual
### Extract IDs of all beggars that are treated as learners
beggar_learners <- as.vector(subset(spongers$id_individual, subset=spongers$Demons_sponging_forage=="no"))
### extract IDs of all beggars treated as demonstrators
beggar_demons <- as.vector(subset(spongers$id_individual, subset=spongers$Demons_sponging_forage=="yes"))
