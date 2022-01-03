# TO DO: Explain what this script is for, and how to use it.
#
# Refer to:
# de_analysis_detailed_look.html
# de_analysis_detailed_look_more.html
#
# TO DO: Add sinteractive command used.
#
library(tools)
library(Matrix)
library(scran)
library(DESeq2)
library(Seurat)
library(MAST)
library(fastTopics)
source("../code/de_analysis_functions.R")

# These variables control the simulations: ns is the number of
# simulations to run; k is the number of topics to simulate; m is the
# number of genes to simulate; num.mc, the length of the Markov chain
# simulated for the DEanalysis; "alpha" is passed to rdirichlet to
# simulate the topic proportions; "outfile" is the file where the
# results of the simulations are stored.
ns      <- 4 # 20
k       <- 2 # 6
m       <- 1e4
num.mc  <- 1000 # 1e4
alpha   <- rep(1,k) # rep(0.01,k)
outfile <- "sims.RData"

# This data structure will be used to store the results of the
# simulations.
res <- vector("list",ns)
names(res) <- paste0("sim",1:ns)

# Repeat for each simulation.
for (i in 1:ns) {
  cat(sprintf("SIMULATION %d\n",i))

  # Simulate counts from a rank-2 Poisson NMF model with parameters
  # chosen to roughly mimic the UMI counts from a single-cell RNA
  # sequencing experiment.
  cat("Generating data set.\n")
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
  cat("Fitting multinomial topic model.\n")
  fit0 <- init_poisson_nmf(X,F = dat$F,L = with(dat,s*L))
  fit <- fit_poisson_nmf(X,fit0 = fit0,numiter = 40,method = "scd",
                         update.loadings = NULL,verbose = "none",
                         control = list(nc = 2))
  fit <- poisson2multinom(fit)

  # Perform a DE analysis *without* shrinking the LFC estimates.
  cat("Performing DE analysis without shrinkage.\n")
  de0 <- de_analysis(fit,X,shrink.method = "none",verbose = FALSE,
                     control = list(ns = num.mc,nc = 2))


  # Perform a second DE analysis using adaptive shrinkage to shrink
  # (and hopefully improve accuracy of) the LFC estimates.
  cat("Performing first DE analysis with shrinkage.\n")
  de1 <- de_analysis(fit,X,shrink.method = "ash",verbose = FALSE,
                     control = list(ns = num.mc,nc = 2))

  # Perform a third DE analysis with the exact same settings as the
  # second, but with a different sequence of pseudorandom numbers.
  # This is intended to verify accuracy of the posterior calculations.
  cat("Performing second DE analysis with shrinkage.\n")
  de2 <- de_analysis(fit,X,shrink.method = "ash",verbose = FALSE,
                     control = list(ns = num.mc,nc = 2))

  # Perform DE analysis using DESeq2. First we prepare the UMI count
  # data for analysis with DESeq2. DESeq2 is called using the settings
  # recommended for single-cell RNA-seq data (see the main DESeq2
  # vignette) To replicate the fastTopics analysis as closely as
  # possible, he LFC estimates are shrunk using ashr.
  cat("Performing DESeq2 analysis.\n")
  Y <- t(X)
  cluster <- factor(apply(dat$L,1,which.max))
  deseq <- DESeqDataSetFromMatrix(Y,data.frame(cluster = cluster),~cluster)
  sizeFactors(deseq) <- calculateSumFactors(Y)
  deseq <- DESeq(deseq,test = "LRT",reduced = ~1,useT = TRUE,minmu = 1e-6,
                 minReplicatesForReplace = Inf)
  deseq <- lfcShrink(deseq,coef = "cluster_2_vs_1",type = "ashr",svalue = TRUE)

  # RUN MAST differential expression analysis via the "FindMarkers"
  # interface in Seurat.
  cat("Performing MAST analysis.\n")
  seurat <- CreateSeuratObject(counts = Y)
  Idents(seurat) <- cluster
  mast <- FindMarkers(seurat,ident.1 = "2",ident.2 = NULL,test.use = "MAST",
                      logfc.threshold = 0,min.pct = 0,verbose = FALSE)
  
  # Store results from the simulations.
  cat("Storing results.\n")
  res[[i]] <- list(data = dat,fit = fit,de0 = de0,de1 = de1,
                   de2 = de2,deseq = deseq,mast = mast)
}

# Write the simulation results to an .RData file.
cat("Saving results to file.\n")
save(list = "res",file = outfile)
resaveRdaFiles(outfile)
