# NOTES:
#
# + Remove ribosomal protein genes?
#
# + Merge topic 1 with topics 3 and 9?
#
library(fastTopics)
library(ggplot2)
library(cowplot)
load("../output/pulseseq/de-pulseseq-merged.RData")
de_merged <- de
rm(de)
p1 <- volcano_plot(de_merged,k = "k7",ymax = 300) # Ciliated cells
p2 <- volcano_plot(de_merged,k = "k1",ymax = 500) # Hillocks
p3 <- volcano_plot(de_merged,k = "k2",ymax = 250) # Tuft, neuroendocrine, ionocytes, proliferating
p4 <- volcano_plot(de_merged,k = "k11",ymax = 500) # Ribosomal protein genes.
p5 <- volcano_plot(de_merged,k = "k3+k9",ymax = 500)
p6 <- volcano_plot(de_merged,k = "k4+k5+k6+k8+k10",ymax = 500)
