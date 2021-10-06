library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load the UMI counts, and the fitted multinomial topic model.
data(pbmc_facs)
counts <- pbmc_facs$counts
fit    <- pbmc_facs$fit

# Perform DE analysis using fitted multinomial topic model.
t0 <- proc.time()
de <- de_analysis(fit,counts)
t1 <- proc.time()
print(t1 - t0)

# For each 
