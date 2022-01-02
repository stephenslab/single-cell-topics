# TO DO: Explain what this script is for, and how to use it.
#
# TO DO: Add sinteractive command used.
#
library(Matrix)
library(fastTopics)
library(tools)
source("../code/de_analysis_functions.R")

# These variables control the simulations: ns is the number of
# simulations to run; k is the number of topics to simulate; m is the
# number of genes to simulate; "outfile" is the file where the results
# of the simulations are stored.
ns      <- 20
k       <- 2
m       <- 10000
outfile <- "sims.RData"

# This data structure will be used to store the results of the
# simulations.
res <- vector("list",ns)
names(res) <- paste0("sim",1:ns)

# Repeat for each simulation.
for (i in 1:ns) {
  cat(sprintf("Running simulation %d:\n",i))

  # Simulate counts from a rank-2 Poisson NMF model with parameters
  # chosen to roughly mimic the UMI counts from a single-cell RNA
  # sequencing experiment.
  cat(" - Generating data set.\n")
  set.seed(i)
  if (k == 2)
    dat <- simulate_twotopic_umi_data(m = m)
  else
    dat <- simulate_manytopic_umi_data(m = m,k = k)
  X   <- dat$X
  F   <- dat$F
  L   <- dat$L

  # Store results from the simulations.
  res[[i]] <- list(data = dat)
}

# Write the simulation results to an RData file.
cat("Saving results to file.\n")
save(list = "res",file = outfile)
resaveRdaFiles(outfile)
