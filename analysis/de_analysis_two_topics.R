library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
source("../code/de_analysis_functions.R")

# Load the results of the simulations with 2 simulated topics.
load("../output/sims/sims-k=2.RData")
res["session.info"] <- NULL

combine_sim_res <- function (res, f)
  do.call("c",lapply(res,function (x) f(x)))

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
  theme_cowplot(font_size = 12)

# Save the plots.
ggsave("../plots/mcmc_scatterplots_sims_k=2.png",
       plot_grid(p1,p2),height = 3,width = 6,dpi = 600)

# Compare K-L divergence measure against fastTopics lfsr, with and
# without adaptive shrinkage.
nonzero_lfc <- combine_sim_res(res,
                 function (x) {
                   i <- which(with(x$fit,F[,2] - F[,1] > -1e-8))
                   return(with(x$dat,abs(F[i,1] - F[i,2]) > 1e-8))
                 })                  
pdat <- data.frame(
  nonzero_lfc = factor(nonzero_lfc),
  noshrink = combine_sim_res(res,
    function (x) {
      i <- which(with(x$fit,F[,2] - F[,1] > -1e-8))
      return(10^(-x$de0$lpval[i,2]))
    }),
  shrink = combine_sim_res(res,
    function (x) {
      i <- which(with(x$fit,F[,2] - F[,1] > -1e-8))
      return(x$de1$svalue[i,2])
     }),
  lfsr = combine_sim_res(res,
    function (x) {
      i <- which(with(x$fit,F[,2] - F[,1] > -1e-8))
      return(x$de1$lfsr[i,2])
     }),
  lkl = combine_sim_res(res,
    function (x) {
      i <- which(with(x$fit,F[,2] - F[,1] > -1e-8))
      return(log10(fastTopics:::min_kl_poisson(x$fit$F)[i,2] + 1e-8))
     }))
p1 <- ggplot(pdat,aes(x = lkl,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  labs(x = "log10 K-L divergence",y = "genes") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = noshrink,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  labs(x = "p-value",y = "genes",title = "without shrinkage") +
  theme_cowplot(font_size = 12)
p3 <- ggplot(pdat,aes(x = shrink,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  labs(x = "s-value",y = "genes",title = "with shrinkage") +
  theme_cowplot(font_size = 12)

# Save the plots.
# TO DO.

# Assess FDR vs. power for all methods.
pdat1 <- create_fdr_vs_power_curve(-pdat$lkl,nonzero_lfc,length.out = 200)
pdat2 <- create_fdr_vs_power_curve(pdat$noshrink,nonzero_lfc,length.out = 200)
pdat3 <- create_fdr_vs_power_curve(pdat$lfsr,nonzero_lfc,length.out = 200)
pdat  <- rbind(cbind(pdat1,method = "kl"),
               cbind(pdat2,method = "noshrink"),
               cbind(pdat3,method = "shrink"))
p <- ggplot(pdat,aes(x = fdr,y = power,color = method)) +
  geom_line(size = 0.65,orientation = "y")  +
  scale_color_manual(values = c("royalblue","firebrick","orange")) +
  theme_cowplot()

# Save the plot.
# TO DO.
