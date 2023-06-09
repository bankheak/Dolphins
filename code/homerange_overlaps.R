# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# HOME RANGE OVERLAPS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

###########################################################################
# PART 1: calculate dyadic home range overlaps ----------------------------

## load all necessary packages
library(tidyverse)
library(raster)
library(viridis)
library(rnaturalearth)
library(ggthemes)
library(patchwork)
library(ggsflabel)
library(ggplot2)
library(rerddapXtracto)
library(sf)
library(sp)
library(rgdal)
library(ggmap)
library(maps)
library(rgdal)
library(spatstat)

# Read in file
sample_data <- read.csv("sample_data.csv")
# Read in map shape file
raster_map <- st_read("Sarasota_Manatee_MPO.shp")
raster_map <- raster(raster_map)

# Extract coordinates
coord_data <- cbind(sample_data[,c(3, 10, 11)]) # Subset Date and Coordinates #
## Format date and year
coord_data$Date <- as.Date(as.character(coord_data$Date), format="%Y-%m-%d")
coord_data$Year <- as.numeric(format(coord_data$Date, format = "%Y"))
## Give descriptive names
colnames(coord_data) <- c("Date", "Latitude", "Longitude", "Year")

# Seperate map per years
years <- unique(coord_data$Year)
coord_years <- list()
for (i in 1:length(years)) {
  coord_years[[i]] <- subset(coord_data, subset=c(coord_data$Year == years[i]))
}    

# Test one year at a time
coord_data <- coord_years[[1]]
## Remove NAs, if any
coord_data = na.omit(coord_data)

# Plot Coordinates on Map
## ggplot theme used for all figures
my_theme <- theme(
  panel.border = element_rect(linewidth =0.5,color="black", fill="transparent"), #white background
  plot.margin=unit(c(2,2,2,2),"mm"), #small margins, first number is top then go clockwise
  panel.background = element_rect(fill = 'white'),
  text=element_text(face="bold", size=8), #general size for all text.looks small here but ok once using ggsave with 600 dpi and default pointsize
  title = element_text(size=rel(1))) # all titles will be a little larger than labels

# North America 
northam <- ne_countries(scale = "medium", returnclass = "sf", continent = "North America")
northam <- subset(northam, admin %in% c("United States of America", "Canada", "Mexico"))

# Florida state only
fl <- maps::map("state", plot = FALSE, fill = TRUE)
fl <- as.data.frame(fl)
fl <- st_as_sf(fl)
fl <- filter(fl, ID == "florida")

# Florida map
FL_map <- ggplotGrob(
  ggplot(northam) +
    geom_sf(fill = "grey95", col = "grey30", lwd = 0.1) +
    geom_sf(data = or, fill = "yellow", col = "grey30", lwd = 0.1) +
    coord_sf(xlim = c(-170, -68), ylim=c(20, 72)) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey95", fill = "transparent", linewidth = 0.5))
  
)

# Now plot
Homerange.plot <- ggplot() +
  annotation_raster(raster_map, xmin = -82.9, xmax = -82.4, ymin = 27, ymax = 27.7) +
  scale_fill_distiller(palette = "Spectral") +
  geom_point(data = coord_data, aes(x = Longitude, y = Latitude)) + 
  coord_quickmap(xlim=c(-82.4, -82.9), ylim=c(27,27.7)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  xlab("Longitude") + ylab("Latitude") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot() +
  geom_tile(data = coord_data) +
  scale_fill_viridis_c(na.value = "transparent", name = "Mean crab fishing effort \n2011-2020") +  
  geom_sf(data = raster_map, col = NA, fill = "grey30") +
  # remove y-axis labels and ticks because they are already shown in panel 1
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(hjust = 0.9)) +
  # add my personal theme
  my_theme &
  # remove x and y titles (because we know it is lat / lon...)
  xlab("") &
  ylab("") &
  # make a few more tweaks to the theme, all panels at the same time
  theme(# move legends to the bottom of each panel
    legend.position = "bottom",
    # make legend title in bold and bigger font
    text=element_text(face="bold", size=8),
    # make legend bars smaller
    legend.key.size = unit(0.8, "lines"),
    # squeeze panels closer together
    panel.spacing = unit(0.5, "lines"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) 
