# Load a few packages.
library(Matrix)
library(DESeq2)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
load("../data/pbmc_purified.RData")

# Load the K = 6 multinomial topic model fit.
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
fit <- poisson2multinom(fit)
