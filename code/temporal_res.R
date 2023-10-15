# 'Multi-network Network-Based Diffusion Analysis

#################################################################################
# TEMPORAL RESOLUTION DEFINITIONS
#################################################################################

# Set working directory here
setwd("../data")

#################################################################################
# PART 1: Divide the data into different resolutions ----------------------------

## load all necessary packages
library(vegan)  
# Run multiple cores for faster computing
library(doParallel)
library(sfsmisc, verbose=F) 

# Read in file and add months
sample_data <- read.csv("sample_data.csv")
list_years <- readRDS("list_years.RData")

# Estimate sampling effort and size for each period
## List years
year_list <- lapply(list_years, function(df) unique(df$Year))
## Get estimate of sampling effort
effort <- as.data.frame(lapply(list_years, function(df) length(unique(df$Date))))
colnames(effort) <- c(1,2)

## Get estimate of population size
unique_ID_year <- as.data.frame(lapply(list_years, function(df) length(unique(df$Code))))
colnames(unique_ID_year) <- c(1,2)

## Get estimate of population size within each HI group
IDbehav_BG <- readRDS("IDbehav_BG.RData")
IDbehav_FG <- readRDS("IDbehav_FG.RData")
IDbehav_SD <- readRDS("IDbehav_SD.RData")

BG <- unique(unlist(sapply(IDbehav_BG, function(df) df$Code[df$HI != 0])))
FG <- unique(unlist(sapply(IDbehav_FG, function(df) df$Code[df$HI != 0])))
SD <- unique(unlist(sapply(IDbehav_SD, function(df) df$Code[df$HI != 0])))

BG_effort <- as.data.frame(lapply(IDbehav_BG, function(df) 
  length(unique(df$Code[df$HI > 0]))))
colnames(BG_effort) <- c(1,2)
FG_effort <- as.data.frame(lapply(IDbehav_FG, function(df) 
  length(unique(df$Code[df$HI > 0]))))
colnames(FG_effort) <- c(1,2)
SD_effort <- as.data.frame(lapply(IDbehav_SD, function(df) 
  length(unique(df$Code[df$HI > 0]))))
colnames(SD_effort) <- c(1,2)

## Compare effort to population size
pop_effort <- as.data.frame(rbind(effort, unique_ID_year, BG_effort, SD_effort, FG_effort)) # Days per year and pop size per year
rownames(pop_effort) <- c('Number of Days Surveyed', 'Number of Individuals', 'Beggars', 'Scavengers/Depredators', 'Fixed Gear Interactors')

# Get all unique Code values in the entire sample_data
all_codes <- unique(sample_data$Code)

# Create a function that counts the IDs in each element
count_instances <- function(df) {
  code_counts <- table(df$Code)
  code_counts <- code_counts[match(all_codes, names(code_counts))] # Add codes to table even if they aren't in that time period
  code_counts[is.na(code_counts)] <- 0 # Replace NAs with 0
  return(code_counts)
}

# Divide resolutions from lowest to highest scale

# -------------------- 22 sets of 1 year increments----------------------------
# Make a list of only 1 year per dataframe
list_years <- split(sample_data, sample_data$Year)
# Apply the count_instances function to each year
instances_per_year <- lapply(list_years, count_instances)
# Convert the list of counts to a data frame
p1y <- do.call(rbind, instances_per_year)
# Transforming into binary matrices
p1y <- as.matrix(p1y); p1y[which(p1y>=1)] = 1; p1y[which(p1y<1)] = 0


# -------------------- 7 sets of 3 year increments----------------------------
# Make a list of 3 years per dataframe
sample_data$ThreeYearIncrement <- cut(sample_data$Year, breaks = seq(min(sample_data$Year), max(sample_data$Year) + 3, by = 3), labels = FALSE)
list_threeyears <- split(sample_data, sample_data$ThreeYearIncrement)
# Apply the count_instances function to each two years
instances_per_threeyear <- lapply(list_threeyears, count_instances)
# Convert the list of counts to a data frame
p3y <- do.call(rbind, instances_per_threeyear)
# Transforming into binary matrices
p3y <- as.matrix(p3y); p3y[which(p3y>=1)] = 1; p3y[which(p3y<1)] = 0


# -------------------- 3 sets of 7 year increments----------------------------
# Make a list of 7 years per dataframe
sample_data$SevenYearIncrement <- cut(sample_data$Year, breaks = seq(min(sample_data$Year), max(sample_data$Year) + 7, by = 7), labels = FALSE)
list_sevenyears <- split(sample_data, sample_data$SevenYearIncrement)
# Apply the count_instances function to each two years
instances_per_sevenyear <- lapply(list_sevenyears, count_instances)
# Convert the list of counts to a data frame
p7y <- do.call(rbind, instances_per_sevenyear)
# Transforming into binary matrices
p7y <- as.matrix(p7y); p7y[which(p7y>=1)] = 1; p7y[which(p7y<1)] = 0


# -------------------- 2 sets of 11 year increments----------------------------
# Make a list of 8 years per dataframe
sample_data$ElevenYearIncrement <- cut(sample_data$Year, breaks = seq(min(sample_data$Year), max(sample_data$Year) + 11, by = 11), labels = FALSE)
list_elevenyears <- split(sample_data, sample_data$ElevenYearIncrement)
# Apply the count_instances function to each two years
instances_per_elevenyear <- lapply(list_elevenyears, count_instances)
# Convert the list of counts to a data frame
p11y <- do.call(rbind, instances_per_elevenyear)
# Transforming into binary matrices
p11y <- as.matrix(p11y); p11y[which(p11y>=1)] = 1; p11y[which(p11y<1)] = 0


###########################################################################
# PART 2: Calculate Whittaker Dissimilarity Index between Time Periods ----------------------------

source("../code/functions.R") # turnover_w function

# Turn over results
t1 = turnover_w(data = p1y, iter = 1000, subseq=F, plot=FALSE)
t3 = turnover_w(data = p3y, iter = 1000, subseq=F, plot=FALSE)
t7 = turnover_w(data = p7y, iter = 1000, subseq=F, plot=FALSE)
t11 = turnover_w(data = p11y, iter = 1000, subseq=F, plot=FALSE)

all = rbind(t1, t3, t7, t11)
all = cbind(c(1, 2, 3, 4), all)

par(mar=c(4,5,4,1))
# Plot the final results. Whisker represent 95%CI generated by the null model. X-axis represent the number of periods and their respective lengths
errbar(x=c(1, 2, 3, 4), y=all[,2], all[,4], all[,5], ylab="Turnover (Averaged Whittaker Dissimilarity)", 
       pch=1, cap=0.02, xaxt='n', xlab="", las=1, cex=1.0, ylim=c(0.30,0.70), xlim=c(1,4), cex.axis=0.8)
axis(1, at=c(1, 2, 3, 4),las=1, labels=c(1, 3, 7, 11), cex.axis=0.7)
mtext(side = 1, "Length of periods (years)", line = 2, font = 1)
axis(3, at=c(1, 2, 3, 4),las=1, labels=c(22, 7, 3, 2), cex.axis=0.7)
mtext(side = 3, "Number of periods", line = 2, font = 1)

# Print final results
all

