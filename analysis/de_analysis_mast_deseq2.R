library(Matrix)
library(Seurat)
library(DESeq2)
library(MAST)
library(fastTopics)
library(ggplot2)
library(cowplot)
source("../code/de_analysis_functions.R")

# Load the results of the simulations when the simulated topic
# proportions are mostly 0 or 1.
load("../output/sims/sims-k=2-alpha=0.01.RData")
res["session.info"] <- NULL

# Most of the simulated cells should be attributed to a single topic.
combine_sim_res <- function (res, f)
  do.call("c",lapply(res,function (x) f(x)))
x <- combine_sim_res(res,function (x) apply(x$dat$L,1,max))
print(mean(x > 0.99))

# Assess accuracy of Monte Carlo estimates (for topic 2 only).
pdat <- data.frame(lfc1 = combine_sim_res(res,function (x) x$de1$postmean[,2]),
                   lfc2 = combine_sim_res(res,function (x) x$de2$postmean[,2]),
                   z1   = combine_sim_res(res,function (x) x$de1$z[,2]),
		   z2   = combine_sim_res(res,function (x) x$de2$z[,2]))
p1 <- ggplot(pdat,aes(x = lfc1,y = lfc2)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "lightsalmon",
              linetype = "dotted") +
  labs(x = "first posterior mean",y = "second posterior mean") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = z1,y = z2)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "lightsalmon",
              linetype = "dotted") +
  labs(x = "first z-score estimate",y = "second z-score estimate") +
  xlim(-45,42) + 
  ylim(-45,42) + 
  theme_cowplot(font_size = 12)

# Save the plots.
ggsave("../plots/mcmc_scatterplots_sims_k=2_alpha=0.01.png",
       plot_grid(p1,p2),height = 3,width = 6,dpi = 600)

# Compare MAST, DESeq2 and fastTopics LFC estimates.
j <- paste0("g",1:10000)
pdat <-
  data.frame(lfc.deseq      = combine_sim_res(res,function (x) x$deseq$log2FoldChange),
             lfc.fasttopics = combine_sim_res(res,function (x) x$de1$postmean[,2]),
             lfc.mast       = combine_sim_res(res,function(x)x$mast[j,"avg_log2FC"]),
             z.deseq        = combine_sim_res(res,
                               function (x) with(x$deseq,log2FoldChange/lfcSE)),
             z.fasttopics   = combine_sim_res(res,function (x) x$de1$z[,2]))
pdat <- transform(pdat,lfc.mast = clamp(lfc.mast,-6,+6))
p1 <- ggplot(pdat,aes(x = lfc.mast,y = lfc.fasttopics)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "lightsalmon",
              linetype = "dotted") +
  labs(x = "MAST",y = "fastTopics",title = "LFC estimates") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = lfc.deseq,y = lfc.fasttopics)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "lightsalmon",
              linetype = "dotted") +
  labs(x = "DESeq2",y = "fastTopics",title = "LFC estimates") +
  theme_cowplot(font_size = 12)
p3 <- ggplot(pdat,aes(x = z.deseq,y = z.fasttopics)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "lightsalmon",
              linetype = "dotted") +
  labs(x = "DESeq2",y = "fastTopics",title = "z-scores") +
  theme_cowplot(font_size = 12)

# Save the plots.
ggsave("../plots/mast_deseq2_sims_scatterplots.png",
       plot_grid(p1,p2,p3,nrow = 1,ncol = 3),
       height = 3,width = 9,dpi = 600)

# Compare MAST, DESeq2 and fastTopics p-values or s-values.
nonzero_lfc <-
  combine_sim_res(res,function (x) with(x$dat,abs(F[,1] - F[,2]) > 1e-8))
pdat <- data.frame(deseq       = combine_sim_res(res,function (x) x$deseq$svalue),
                   fasttopics  = combine_sim_res(res,function (x) x$de1$svalue[,2]),
                   mast        = combine_sim_res(res,function (x) x$mast[j,"p_val"]),
                   nonzero_lfc = factor(nonzero_lfc))
p4 <- ggplot(pdat,aes(x = mast,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  labs(x = "p-value",y = "genes",title = "MAST") +
  theme_cowplot(font_size = 12)
p5 <- ggplot(pdat,aes(x = deseq,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  labs(x = "s-value",y = "genes",title = "DESeq2") +
  theme_cowplot(font_size = 12)
p6 <- ggplot(pdat,aes(x = fasttopics,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  labs(x = "s-value",y = "genes",title = "fastTopics") +
  theme_cowplot(font_size = 12)

# Save the plots.
ggsave("../plots/mast_deseq2_sims_pvalues.eps",
       plot_grid(p4,p5,p6,nrow = 1,ncol = 3),
       height = 2.75,width = 9)

# Assess FDR vs. power for all methods.
pdat1 <- create_fdr_vs_power_curve(pdat$deseq,nonzero_lfc,length.out = 200)
pdat2 <- create_fdr_vs_power_curve(combine_sim_res(res,function (x) x$de1$lfsr[,2]),
                                   nonzero_lfc,length.out = 200)
pdat3 <- create_fdr_vs_power_curve(pdat$mast,nonzero_lfc,length.out = 200)
pdat  <- rbind(cbind(pdat1,method = "deseq2"),
               cbind(pdat2,method = "fastTopics"),
               cbind(pdat3,method = "mast"))
p <- ggplot(pdat,aes(x = fdr,y = power,color = method)) +
  geom_line(size = 0.65,orientation = "y")  +
  scale_color_manual(values = c("dodgerblue","darkorange","darkblue")) +
  theme_cowplot(font_size = 12)

# Save the plot.
ggsave("../plots/mast_deseq2_sims_fdr_vs_power.eps",p,height = 3,width = 4)
