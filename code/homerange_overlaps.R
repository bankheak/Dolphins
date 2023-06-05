# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# HOME RANGE OVERLAPS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: calculate dyadic home range overlaps ----------------------------

## load all necessary packages
library(raster)
library(ggplot2)
library(rerddapXtracto)
library(sf)

# Read in file
sample_data <- read.csv("sample_data.csv")
# Read in map shape file
land_sf <- st_read("Sarasota_Manatee_MPO.shp")

# Extract coordinates
coord_data <- cbind(sample_data[,c(3, 10, 11)]) # Subset Date and Coordinates #
## Format date and year
coord_data$Date <- as.Date(as.character(coord_data$Date), format="%Y-%m-%d")
coord_data$Year <- as.numeric(format(coord_data$Date, format = "%Y"))
## Give descriptive names
colnames(coord_data) <- c("Date", "Latitude", "Longitude", "Year")

# Seperate map per years
years <- unique(coord_data$Year)
list_years <- list()
for (i in 1:length(years)) {
  list_years[[i]] <- subset(coord_data, subset=c(coord_data$Year == years[i]))
}    

# Test one year at a time
coord_data <- list_years[[1]]

# Plot Coordinates on Map
# now plot
Homerange.plot <- ggplot() +
  # plot sighting point locations 
  geom_tile(data = land_sf, aes(x = Longitude, y = Latitude)) +
  scale_fill_distiller(palette = "Spectral") +
  geom_point(data = coord_data, aes(x = Longitude, y = Latitude)) + 
  coord_quickmap(xlim=c(-82, -83), ylim=c(27,28)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  xlab("Longitude") + ylab("Latitude") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot() +
geom_sf(data = land_sf, col = "grey10", fill = "grey15") +
 #scale_color_viridis_c(guide = F) +
  #mon_theme +
  geom_point(data = coord_data, aes(x = Longitude, y = Latitude)) + 
  coord_sf(expand = T) +
  #coord_quickmap(xlim=c(-82, -83), ylim=c(27,28)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(axis.text.y  = element_text(angle = 90, hjust = 0.5)) 


