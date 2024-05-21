# 'Multi-network Network-Based Diffusion Analysis

#########################################################################################
# HOME RANGE OVERLAPS
#########################################################################################

# Set working directory here
setwd("../data")

#########################################################################################
# PART 1: calculate dyadic home range overlaps of individuals----------------------------

## load all necessary packages
library(MASS)
library(sf) # Convert degrees to meters
library(sp) # Creates a SpatialPointsDataFrame by defining the coordinates
library(adehabitatHR) # Caluculate MCPs and Kernel density 
library(scales) # Helps make polygons partly transparent using the alpha argument
library(ggmap) # Download tiles using ggmap
library(viridis) # Color pallette
library(gridExtra) # grid.arrange function
library(ggplot2)
library(patchwork) # align different plots
library(dplyr)
library(adehabitatHR)
library(ggplot2)
library(magrittr)
library(leaflet)

# Read in file
list_years <- readRDS("list_years.RData")

# Aggregate list into one homerange overlap matrix
list_years_df <- merge(list_years[[1]], list_years[[2]], all = T)

# Transform coordinate data into a Spatial Points Dataframe in km
create_coord_data <- function(df) {
  
    # Extract IDs and coordinates
    ids <- df$Code
    coordinates <- df[, c("StartLon", "StartLat")]
    
    # Create a SpatialPointsDataFrame with coordinates
    coords_sp <- SpatialPointsDataFrame(coords = coordinates, data = data.frame(id = ids))
    
    # Set CRS to WGS84
    proj4string(coords_sp) <- CRS("+proj=longlat +datum=WGS84")
    
    # Transform to a UTM CRS that uses km as the unit
    coords_sp_utm <- spTransform(coords_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
    
    return(coords_sp_utm)
}

dolph.sp <- create_coord_data(list_years_df)

# Visualize data extent
vis.sf <- function(dolph.sp) {
dolp.sf <- st_as_sf(dolph.sp)
  return(dolp.sf)
}

dolph.sf <- vis.sf(dolph.sp)

vis_coord <- lapply(seq_along(dolph.sf), function (i) {ggplot(dolph.sf[[i]]) +
                      geom_sf(aes(color = "Dolphins"), size = 2, alpha = 0.5) +
                      theme_bw() +
                      labs(title = "Distribution") +
                      scale_color_manual(values = c("Dolphins" = "blue"))})
wrap_plots(vis_coord, nrow = 1)

# Look into what bandwidth 
value <- 1
period <- 1
n <- length(dolph.sp[[period]]@coords[, value]) 
bw_scott <- (4 / (3 * n))^(1/5) * sd(dolph.sp[[period]]@coords[, value])  # Scott's rule
## Use the calculated bandwidth in density estimation
density_estimate <- density(dolph.sp[[period]]@coords[, value], bw = bw_scott)
## Plot the density estimate
plot(density_estimate)

# Calculate the maximum distance between points
max(sqrt((dolph.sp[[period]]@coords[,1] - min(dolph.sp[[period]]@coords[,1]))^2 + 
           (dolph.sp[[period]]@coords[,2] - min(dolph.sp[[period]]@coords[,2]))^2))

# Use the calculated extent in kernelUD
kernel <- kernelUD(dolph.sp, h = 1000)

# Get centroid points for individuals
kernel <- readRDS("kernel.RData")

###########################################################################
# PART 2: Calculate Dyadic HRO Matrix: HRO = (Rij/Ri) * (Rij/Rj)------------------------------------------------------------

# Calculate kernel overlap values
kov <- kerneloverlaphr(kernel, method = "HR", lev = 95)

# Save HRO
saveRDS(kov, "kov.RDS")

###########################################################################
# PART 3: Plot HRO for HI Dolphins ------------------------------------------------------------

# Make coordinate data
dolph.data <- lapply(list_years, function (df) {
  # Read in only ID, x, and y
  df <- df[, c("Code", "StartLon", "StartLat")]
  # Create a SpatialPointsDataFrame by defining the coordinates
  coordinates(df) <- c("StartLon", "StartLat")
  # Set the coordinate reference system (CRS) 
  # More information on CRS here: 
  # https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf 
  proj4string(df) <- CRS( "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs" )
  
  df
  })

write.csv(dolph.data[[1]], "dolph.data.csv")

# Calculate MCPs for each dolphin
dolph.mcp <- lapply(dolph.data, function (df) mcp(df, percent = 100))
plot(dolph.data[[1]], col = as.factor(dolph.data[[1]]@data$id), pch = 16) 
plot(dolph.mcp[[1]], col = alpha(1:5, 0.5), add = TRUE)

# Transform the point and MCP objects. 
dolph.spgeo <- lapply(dolph.data, function (df) spTransform(df, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")))
dolph.mcpgeo <- lapply(dolph.mcp, function (df) spTransform(df, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")))

# Google tiles (requires a key first)
## Guide on getting a key: https://www.r-bloggers.com/geocoding-with-ggmap-and-the-google-api/
register_google(key = "AIzaSyAgFfxIJmkL8LAWE7kHCqSqKBQDvqa9umI")
mybasemap <- get_map(location = c(lon = mean(dolph.spgeo[[1]]@coords[,1]),
                                  lat = mean(dolph.spgeo[[1]]@coords[,2])),
                     source = "google",
                     zoom = 10,
                    maptype = 'satellite')

# Turn the spatial data frame of points into just a dataframe for plotting in ggmap
dolph.geo <- data.frame(dolph.spgeo[[1]]@coords, 
                          id = dolph.spgeo[[1]]@data$Code )

mymap.hr <- ggmap(mybasemap) + 
  geom_polygon(data = fortify(dolph.mcpgeo[[1]]),  
               # Polygon layer needs to be "fortified" to add geometry to the dataframe
               aes(long, lat, colour = id, fill = id),
               alpha = 0.3) + # alpha sets the transparency
  geom_point(data = dolph.geo, 
             aes(x = coords.x1, y = coords.x2, colour = id))  +
  labs(x = "Longitude", y = "Latitude") + theme(legend.position = "none")
mymap.hr

# Find HI events among individuals
ID_HI <- lapply(list_years, function(df) {
  subset_df <- subset(df, subset = c(df$ConfHI != 0))
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
dolph.sp_HI <- lapply(ID_HI, function (df) create_coord_data(df))

# Calculate MCPs for each HI dolphin
HI.mcp <- lapply(dolph.sp_HI, function (df) mcp(df, percent = 95))

# Transform the point and MCP objects. 
HI.spgeo <- lapply(ID_HI, function (df) spTransform(df, CRS("+proj=longlat")))
HI.mcpgeo <- lapply(HI.mcp, function (df) spTransform(df, CRS("+proj=longlat")))

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

