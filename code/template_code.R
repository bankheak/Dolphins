# 'Multi-network Network-Based Diffusion Analysis

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins")


# PART 1: calculate dyadic home range overlaps ----------------------------

## load all necessary packages
require(rmarkdown)  # “Knit” button (Ctrl+Shift+K) displays preview
require(ggplot2)  # Plot designs



# PART 2: applying NBDA to 'begging' data ---------------------------------

## load all necessary packages
install.packages("devtools")
require(devtools)
install_github("whoppitt/NBDA")
require(NBDA)

## Read in all network data and turn into matrix

# Vertical network
vert.data <- read.csv("../data/social_vertical.csv", row.names=1, header=TRUE)
vert.data <- as.matrix(vert.data)

# Horizontal network
horz.data <- read.csv("../data/social_horizontal.csv", row.names=1, header=TRUE)
horz.data <- as.matrix(horz.data)

# Ecological network
ecol.data <- read.csv("../data/HR_overlaps.csv", row.names=1, header=TRUE)
ecol.data <- as.matrix(ecol.data)

# Relatedness network
relate.data <- read.csv("../data/relatedness.csv", row.names=1, header=TRUE)
relate.data <- as.matrix(relate.data)

# read ILVs
ILV.data <- read.csv("../data/ILVs.csv", header=TRUE, sep=",")
ILV.data[ILV.data=="<NA>"]=NA

