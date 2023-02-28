library(tidyverse)
library(ggplot2)
library(readr)
library(maps)
library(viridis)
library(rcartocolor)

#input data
mapdata = read.csv("taraocean_count.csv", header=T, row.names = 1)
world <- map_data("world")

#plot map with coordinate
ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey") +
  geom_point(data=mapdata, aes(x=long, y=lat, size=cpm1, color=cpm1),stroke=F, alpha=0.8)+
  scale_size_continuous(range=c(0.1,7))+
  scale_color_carto_c(palette = "BurgYl")+
  guides( colour = guide_legend()) +
  theme(
    legend.position = "bottom",
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#ffffff", color = NA), 
    panel.background = element_rect(fill = "#ffffff", color = NA), 
    legend.background = element_rect(fill = "#ffffff", color = NA)
  )

#plot map without coordinate
ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey") +
  geom_point(data=mapdata, aes(x=long, y=lat, size=cpm1, color=cpm1),stroke=F, alpha=0.8)+
  scale_size_continuous(range=c(0.1,7))+
  scale_color_carto_c(palette = "BurgYl")+ #Earth BurgYl
    theme_void() + 
  guides( colour = guide_legend()) +
  theme(
    legend.position = "bottom",
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#ffffff", color = NA), 
    panel.background = element_rect(fill = "#ffffff", color = NA), 
    legend.background = element_rect(fill = "#ffffff", color = NA)
  )
