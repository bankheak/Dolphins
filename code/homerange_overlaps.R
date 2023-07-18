# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# HOME RANGE OVERLAPS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: calculate dyadic home range overlaps ----------------------------

## load all necessary packages
library(sf) # Convert degrees to meters
library(sp) # Creates a SpatialPointsDataFrame by defining the coordinates
library(adehabitatHR) # Caluculate MCPs
library(scales) # Helps make polygons partly transparent using the alpha argument
library(ggmap) # Download tiles using ggmap
library(viridis) # Color pallette
library(gridExtra) # grid.arrange function
library(ggplot2)
library(adehabitatHR) # Kernel density 
library(rgdal) # Overlap

# Read in file
sample_data <- read.csv("sample_data.csv")

# Create background map using ggmap
mybasemap <- ggmap(get_stamenmap(bbox = c(left = -83, bottom = 27, right = -82, top = 28)))

# Extract coordinates
coord_data <- cbind(sample_data[,c('Date', 'StartLat', 'StartLon', 'Code', 'subYear')]) # Subset Date and Coordinates #
## Format date and year
coord_data$Date <- as.Date(as.character(coord_data$Date), format="%Y-%m-%d")
## Give descriptive names
colnames(coord_data) <- c("date", "y", "x", "id", "subyear")

# Seperate map per years
years <- unique(coord_data$subyear)
coord_years <- list()
for (i in 1:length(years)) {
  coord_years[[i]] <- subset(coord_data, subset=c(coord_data$subyear == years[i]))
}    

# Test one year at a time
coord_data <- coord_years[[1]]

# Eliminate IDs with less than 5 locations
coord_data <- subset(coord_data, subset=c(coord_data$id != "None"))
coord_data <- coord_data[!is.na(coord_data$x) & !is.na(coord_data$y),]

ID <- unique(coord_data$id)
obs_vect <- NULL
for (i in 1:length(ID)) {
  obs_vect[i]<- sum(coord_data$id == ID[i])
}
sub <- data.frame(ID, obs_vect)
sub <- subset(sub, subset=c(sub$obs_vect > 10))
coord_data <- subset(coord_data, coord_data$id %in% c(sub$ID))

# Only include three columns (id, x, and y coordinates) for making MCP's
dolph.sp <- coord_data[, c("id", "y", "x")] 

# Create a simple feature data frame (sf)
coord_data_sf <- st_as_sf(dolph.sp, coords = c("x", "y"), crs = 4326)

# UTM zone for study area
dolph.sf <- st_transform(coord_data_sf, crs = paste0("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
# Extract coordinates (latitude and longitude) and create new columns
dolph.sp$x <- st_coordinates(dolph.sf)[, 1]
dolph.sp$y <- st_coordinates(dolph.sf)[, 2]
# Remove two rows with NA's
dolph.sp <- dolph.sp[!is.na(dolph.sp$x) & !is.na(dolph.sp$y),]

coordinates(dolph.sp) <- c("x", "y")

# Set the initial CRS for data to WGS84 (latitude and longitude)
proj4string(dolph.sp) <- CRS( "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs" )

# Create kernel estimates for each id
kernel.lscv <- kernelUD(dolph.sp, h = "LSCV")  # LSCV = least squares cross validation

# Get the names or IDs of the individuals that are of concern
selected_individuals <- c("F101", "F191", "F192", "FB15", "FB35", "FB41", "FB93")

# Repeat code above to calculate appropriate bandwidth for IDs of concern
concern.sp <- coord_data[, c("id", "y", "x")] 
concern.sp <- subset(concern.sp, id %in% selected_individuals) # Only include selected ids

concern_sf <- st_as_sf(concern.sp, coords = c("x", "y"), crs = 4326)
concern.sf <- st_transform(concern_sf, crs = paste0("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
concern.sp$x <- st_coordinates(concern.sf)[, 1]
concern.sp$y <- st_coordinates(concern.sf)[, 2]
concern.sp <- concern.sp[!is.na(concern.sp$x) & !is.na(concern.sp$y),]

coordinates(concern.sp) <- c("x", "y")

proj4string(concern.sp) <- CRS( "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs" )

# Now recalculate kernel estimates for each id
kernel.con <- kernelUD(concern.sp, h = "LSCV")
plotLSCV(kernel.con) # it looks like a bandwidth of 1000 will be good enough

kernel.lscv <- kernelUD(dolph.sp, h = 1000)

# 95% of estimated distribution
dolph.kernel.poly <- getverticeshr(kernel.lscv, percent = 95)  
print(dolph.kernel.poly)  # returns the area of each polygon

###########################################################################
# PART 2: Plot HRO ------------------------------------------------------------

# Calculate MCPs for each dolphin
dolph.mcp <- mcp(dolph.sp, percent = 95)

# Convert dolph.sp and dolph.mcp to sf objects
dolph.sp_sf <- st_as_sf(dolph.sp)
dolph.mcp_sf <- st_as_sf(dolph.mcp)

# Convert ggmap basemap to raster
mybasemap_raster <- as.raster(mybasemap)

# Plot the ggmap basemap as a raster and overlay the two sf objects
ggplot() +
  annotation_raster(mybasemap_raster, xmin = -83, xmax = -82, ymin = 27, ymax = 28) +
  geom_sf(data = dolph.mcp_sf, aes(fill = as.factor(id)), alpha = 0.5) +
  geom_sf(data = dolph.sp_sf, aes(color = factor(id)), size = 3, shape = 16) +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8, direction = -1) +  # Use viridis color palette
  scale_color_manual(values = rainbow(length(unique(dolph.sp_sf$id)))) +
  theme_minimal() +
  labs(fill = "MCP", color = "Dolphin number") +
  coord_sf() + theme(legend.position = "none")  # Remove the legend

###########################################################################
# PART 2: Calculate Dyadic HRO Matrix: HRO = (Rij/Ri) * (Rij/Rj)------------------------------------------------------------

# Get HRO 
kov <- kerneloverlaphr(kernel.lscv, method="HR", lev=95)
