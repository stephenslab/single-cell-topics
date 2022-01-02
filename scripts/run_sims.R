# TO DO: Explain what this script is for, and how to use it.
#
# Refer to:
# de_analysis_detailed_look.html
# de_analysis_detailed_look_more.html
#
# TO DO: Add sinteractive command used.
#
library(Matrix)
library(fastTopics)
library(tools)
source("../code/de_analysis_functions.R")

# These variables control the simulations: ns is the number of
# simulations to run; k is the number of topics to simulate; m is the
# number of genes to simulate; "alpha" is passed to rdirichlet to
# simulate the topic proportions; "outfile" is the file where the
# results of the simulations are stored.
ns      <- 20
k       <- 2 # 6
m       <- 10000
alpha   <- rep(1,k) # rep(0.01,k)
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
    dat <- simulate_twotopic_umi_data(m,alpha = alpha)
  else
    dat <- simulate_manytopic_umi_data(m,k = k,alpha = alpha)
  X <- dat$X

  # Fit a multinomial topic model to the simulated UMI count data. To
  # simplify evaluation, L is assumed to be known, and fix them to
  # their ground-truth values. In this way, the only error that can
  # arise is in the estimates of F.
  cat(" - Fitting multinomial topic model.\n")
  fit0 <- init_poisson_nmf(X,F = dat$F,L = with(dat,s*L))
  fit <- suppressMessages(
           fit_poisson_nmf(X,fit0 = fit0,numiter = 40,method = "scd",
                           update.loadings = NULL,verbose = "none",
                           control = list(nc = 2)))
  fit <- poisson2multinom(fit)

  # Store results from the simulations.
  res[[i]] <- list(data = dat,fit = fit)

  stop()
}

# Write the simulation results to an RData file.
cat("Saving results to file.\n")
save(list = "res",file = outfile)
resaveRdaFiles(outfile)
