# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# HOME RANGE OVERLAPS
###########################################################################

# Set working directory here
setwd("../data")

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
list_sexage_years <- readRDS("list_sexage_years.RData")

# Make a function that calculates all of the following code with and without sex and age
create_coord_data <- function(list_years) {

# Make a list of years
coord_data_list <- list_years

# Process the coord_data_list
dolph.sp <- lapply(coord_data_list, function(df) {
  # Extract IDs and coordinates
  ids <- df$Code
  coordinates <- df[, c("StartLon", "StartLat")]
  
  # Convert to data frame
  ids_df <- data.frame(id = ids)
  
  # Create a SpatialPointsDataFrame with coordinates
  coords_sp <- SpatialPointsDataFrame(coords = coordinates, data = ids_df)
  
  # Set CRS and transform to UTM
  proj4string(coords_sp) <- CRS("+proj=longlat +datum=WGS84")
  coords_sp_utm <- spTransform(coords_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
})
return(dolph.sp)
}

dolph.sp <- create_coord_data(list_years)
dolph.sp_sexage <- create_coord_data(list_sexage_years)

# Visualize data extent
vis.sf <- function(dolph.sp) {
dolph.sf <- lapply(dolph.sp, function (df) {st_as_sf(df)})
ggplot(dolph.sf[[3]]) +
  geom_sf(aes(color = "Data Points"), size = 2, alpha = 0.5) +
  theme_bw() +
  labs(title = "Distribution of Data Points") +
  scale_color_manual(values = c("Data Points" = "blue"))
return(dolph.sf)
}

dolph.sf <- vis.sf(dolph.sp)
dolph.sf_sexage <- vis.sf(dolph.sp_sexage)

# Calculate kernel values
create_kernel <- function(dolph.sp) {
kernel <- lapply(dolph.sp, function(sp_obj) {
  kernelUD(sp_obj, h = 10000)
})
return(kernel)
}

kernel <- create_kernel(dolph.sp)
kernel_sexage <- create_kernel(dolph.sp_sexage)

# Create area of each polygon
year <- 5
dolph.kernel.poly <- getverticeshr(kernel[[year]], percent = 95)
print(dolph.kernel.poly)  

###########################################################################
# PART 2: Calculate Dyadic HRO Matrix: HRO = (Rij/Ri) * (Rij/Rj)------------------------------------------------------------

# Calculate kernel overlap values
create_kov <- function(kernel) {
kov <- lapply(kernel, function(kern) {
  kerneloverlaphr(kern, method = "HR", lev = 95)
})
return(kov)
}

kov <- create_kov(kernel)
kov_sexage <- create_kov(kernel_sexage)

saveRDS(kov, "kov.RDS")
saveRDS(kov_sexage, "kov_sexage.RDS")

###########################################################################
# PART 3: Plot HRO for HI Dolphins ------------------------------------------------------------

# Find HI events among individuals
ID_HI <- lapply(coord_data_list, function(df) {
  subset_df <- subset(df, subset = df$ConfHI != 0)
  subset_df <- subset_df[, c('StartLat', 'StartLon', 'Code')]
  
  # Make sure there are at least 5 relocations
  ID_df <- unique(subset_df$Code)
  
  obs_vect <- numeric(length(ID_df))
  for (i in seq_along(ID_df)) {
    obs_vect[i] <- sum(subset_df$Code == ID_df[i])
  }
  
  sub <- data.frame(ID_df, obs_vect)
  sub <- subset(sub, subset = obs_vect > 4)
  
  subset_df <- subset_df[subset_df$Code %in% sub$ID_df, ]
  
  return(subset_df)
})

# Recalculate Coordinate data
HI.sp <- lapply(ID_HI, function(df) {
  # Extract IDs and coordinates
  ids <- df$Code
  coordinates <- df[, c("StartLon", "StartLat")]
  
  # Convert to data frame
  ids_df <- data.frame(id = ids)
  
  # Create a SpatialPointsDataFrame with coordinates
  coords_sp <- SpatialPointsDataFrame(coords = coordinates, data = ids_df)
  
  # Set CRS and transform to UTM
  proj4string(coords_sp) <- CRS("+proj=longlat +datum=WGS84")
  coords_sp_utm <- spTransform(coords_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  
  # Calculate kernel density estimates
  kernel_obj <- kernelUD(coords_sp_utm, h = 1000)
  
  # Add the kernel density estimates to the SpatialPointsDataFrame
  coords_sp_utm$estUD <- kernel_obj$estUD
  
  coords_sp_utm
})

# Kernel estimate
HI.kern <- lapply(HI.sp, function(sp_obj) {
  kernelUD(sp_obj, h = 500)
})

HI.kernel.poly <- getverticeshr(HI.kern[[year]], percent = 95)

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
