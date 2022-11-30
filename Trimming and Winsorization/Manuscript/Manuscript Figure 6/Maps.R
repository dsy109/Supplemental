library(maps)
library(urbnmapr)
library(ggplot2)
library(tidyverse)
setwd("/Users/kcheng/Desktop/Real Data/Map Plot")
county.data <- read_csv("County.csv")
###############################################################
###############################################################
### Combine County.data with Longitude and Latitude ###
farm_data <- left_join(county.data, counties, by = "county_fips") 
### Draw Maps ###
farm_data %>%
  ggplot(aes(long, lat, group = group, fill = AdjCount)) +
  geom_polygon(color = "azure4" , size=0.15) +
  scale_fill_gradientn(colors=rev(terrain.colors(5))) +
  geom_polygon(data = states, mapping = aes(long, lat, group = group),
               fill = NA, color = "black") +
  theme(legend.title = element_text(size = 25 ,face = "bold"),
        legend.position = c(0.95,0.18),
        legend.text = element_text(size=25,face = "bold"),
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(fill = "Count")
  
