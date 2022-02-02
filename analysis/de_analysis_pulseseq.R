library(fastTopics)
library(ggplot2)
library(cowplot)
load("../output/pulseseq/de-pulseseq-merged.RData")
de_merged <- de
rm(de)
p3 <- volcano_plot(de_merged,k = "k7",ymax = 500)
