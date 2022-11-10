# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# DYATIC SOCIAL ASSOCIATIONS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/OneDrive/Documents/Homework/OSU/R Homework")

###########################################################################
# PART 1: Gambit of the group ---------------------------------------------

# Load all necessary packages
require(dplyr)

# Read file in
group_data <- read.csv("group.data.csv")

# Make date into a date class
group_data$Date <- as.Date(as.character(group_data$Date), format="%m/%d/%Y")
group_data <- rbind(group_data, group_data[1:5,])  # Add data for better use

# Group each individual by date and group
check <-group_data |>
  group_by(Individual, Date, Group_number) |>
  summarize(indi_group_date = n())

# Create association index


# Calculate simple-ratio index
