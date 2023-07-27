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
library(adehabitatHR) # Caluculate MCPs and Kernel density 
library(scales) # Helps make polygons partly transparent using the alpha argument
library(ggmap) # Download tiles using ggmap
library(viridis) # Color pallette
library(gridExtra) # grid.arrange function
library(ggplot2)
library(rgdal) # Overlap

# Read in file
sample_data <- read.csv("sample_data.csv")
list_years <- readRDS("list_years.RData")

## Test one year at a time
year <- 7
coord_data <- list_years[[year]]

# Extract coordinates
coord_data <- cbind(coord_data[,c('Date', 'StartLat', 'StartLon', 'Code', 'Year', 'ConfHI')]) # Subset Date and Coordinates #
## Format date and year
coord_data$Date <- as.Date(as.character(coord_data$Date), format="%Y-%m-%d")
## Give descriptive names
colnames(coord_data) <- c("date", "y", "x", "id", "year", "HI")

# Only include three columns (id, x, and y coordinates) for making MCP's
dolph.sp <- coord_data[, c("id", "y", "x")] 

# Create a simple feature data frame (sf)
coord_data_sf <- st_as_sf(dolph.sp, coords = c("x", "y"), crs = 4326)

# UTM zone for study area
dolph.sf <- st_transform(coord_data_sf, crs = paste0("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))

# Extract coordinates (latitude and longitude) and create new columns
dolph.sp$x <- st_coordinates(dolph.sf)[, 1]
dolph.sp$y <- st_coordinates(dolph.sf)[, 2]

coordinates(dolph.sp) <- c("x", "y")

# Set the initial CRS for data to WGS84 (latitude and longitude)
proj4string(dolph.sp) <- CRS( "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs" )

# Create kernel estimates for each id
kernel.lscv <- kernelUD(dolph.sp, h = "LSCV")  # LSCV = least squares cross validation

# Get the names or IDs of the individuals that are of concern
selected_individuals <- c("BEGR", "1091", "1316", "BRN1")

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

# 95% of estimated distribution
kernel <- kernelUD(dolph.sp, h = 500)
dolph.kernel.poly <- getverticeshr(kernel, percent = 95)  
print(dolph.kernel.poly)  # returns the area of each polygon

###########################################################################
# PART 2: Plot HRO for HI Dolphins ------------------------------------------------------------

# Find HI events among individuals
ID_HI <- subset(coord_data, subset=c(coord_data$HI != 0))
ID_HI <- ID_HI[,c('y', 'x', 'id')] 

# Make sure there are at least 5 relocations
ID <- unique(ID_HI$id)
obs_vect <- NULL
for (i in 1:length(ID)) {
  obs_vect[i]<- sum(ID_HI$id == ID[i])
}
sub <- data.frame(ID, obs_vect)
sub <- subset(sub, subset=c(sub$obs_vect > 4))
ID_HI <- subset(ID_HI, ID_HI$id %in% sub$ID)

# Recalculate Coordinate data
ID_HI_sf <- st_as_sf(ID_HI, coords = c("x", "y"), crs = 4326)
HI.sf <- st_transform(ID_HI_sf, crs = paste0("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
ID_HI$x <- st_coordinates(HI.sf)[, 1]
ID_HI$y <- st_coordinates(HI.sf)[, 2]

ID_HI <- ID_HI[!is.na(ID_HI$x) & !is.na(ID_HI$y),]

coordinates(ID_HI) <- c("x", "y")

proj4string(ID_HI) <- CRS( "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs" )

# Kernel estimate
HI.kern <- kernelUD(ID_HI, h = 500)
HI.kernel.poly <- getverticeshr(HI.kern, percent = 95)

# Plot kernel density
colors <- rainbow(length(unique(ID_HI$id)))
individuals <- unique(HI.kernel.poly@data$id)
## Match each individual to a color
individual_color <- colors[match(individuals, unique(HI.kernel.poly@data$id))]
## Match the color for each home range polygon
color <- individual_color[match(HI.kernel.poly@data$id, individuals)]
## Plot the home range polygons with colors
plot(HI.kernel.poly, col = color)

# Calculate MCPs for each HI dolphin
HI.mcp <- mcp(ID_HI, percent = 95)

# Transform the point and MCP objects. 
HI.spgeo <- spTransform(ID_HI, CRS("+proj=longlat"))
HI.mcpgeo <- spTransform(HI.mcp, CRS("+proj=longlat"))

# Turn the spatial data frame of points into just a dataframe for plotting in ggmap
HI.geo <- data.frame(HI.spgeo@coords, 
                          id = HI.spgeo@data$id )

# Create background map using ggmap
mybasemap <- get_stamenmap(bbox = c(left = -83, bottom = 27, right = -82, top = 28))

# Plot HI ids
mymap.hr <- ggmap(mybasemap) + 
  geom_polygon(data = fortify(HI.mcpgeo),  
               # Polygon layer needs to be "fortified" to add geometry to the dataframe
               aes(long, lat, colour = id, fill = id),
               alpha = 0.3) + # alpha sets the transparency
  geom_point(data = HI.geo, 
             aes(x = x, y = y, colour = id))  +
  theme(legend.position = c(0.15, 0.80)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_manual(name = "Dolphin ID", 
                    values = colors,
                    breaks = individuals) +
  scale_colour_manual(name = "Dolphin ID", 
                      values = colors,
                      breaks = individuals)

# Categorize ConfHI to IDs
ID_HI <- subset(coord_data, subset=c(coord_data$HI != 0))

HI_type <- as.matrix(table(ID_HI$id, ID_HI$HI))
HI_type <- data.frame(HI_type)
colnames(HI_type) <- c("Code", "ConfHI", "Freq")

HI <- subset(HI_type, subset=c(HI_type$Freq != 0))
HI <- subset(HI, HI$Code %in% individuals)

###########################################################################
# PART 2: Calculate Dyadic HRO Matrix: HRO = (Rij/Ri) * (Rij/Rj)------------------------------------------------------------

# Get HRO 
kernel <- kernelUD(dolph.sp, h = 1000)

kov <- kerneloverlaphr(kernel, method="HR", lev=95)

saveRDS(kov, "kov.RDS")
