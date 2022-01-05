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
