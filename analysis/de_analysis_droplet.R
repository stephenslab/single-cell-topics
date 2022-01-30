# Load the packages and used in the analysis.
library(fastTopics)
library(ggplot2)
library(cowplot)

# Load the results of the differential expression analysis.
load("../output/droplet/de-droplet.RData")

# Summarize the results of the DE analysis in volcano plots.
p1 <- volcano_plot(de_merged,ymax = 100,k = "k1") # Goblet cells
p2 <- volcano_plot(de_merged,ymax = 100,k = "k2") # Basal cells
p4 <- volcano_plot(de_merged,ymax = 200,k = "k4") # Hillocks
p5 <- volcano_plot(de_merged,ymax = 200,k = "k5+k7") # Club cells
p6 <- volcano_plot(de_merged,ymax = 100,k = "k6") # Ciliated cells
p7 <- volcano_plot(de,k = "k5",ymax = 200) # Other club cells
p8 <- volcano_plot(de,k = "k7",ymax = 200) # Scgb1a1+ club cells
plot_grid(p7,p8)
