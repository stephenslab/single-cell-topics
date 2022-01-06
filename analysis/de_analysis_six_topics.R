library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
source("../code/de_analysis_functions.R")
set.seed(1)

# Load the results of the simulations.
load("../output/sims/sims-k=6.RData")

# Before comparing the methods, we first assess accuracy of the Monte
# Carlo computations by comparing the estimates from two independent
# MCMC runs.
res["session.info"] <- NULL
pdat <- data.frame(lfc1 = combine_sim_res(res,function (x) c(x$de1$postmean)),
                   lfc2 = combine_sim_res(res,function (x) c(x$de2$postmean)),
                   z1   = combine_sim_res(res,function (x) c(x$de1$z)),
	           z2   = combine_sim_res(res,function (x) c(x$de2$z)))
i <- which(abs(pdat$z1) <= 1 & abs(pdat$z2) <= 1)
j <- which(abs(pdat$z1) > 1 | abs(pdat$z2) > 1)
i <- sample(i,1e5)
pdat <- pdat[c(i,j),]
p1 <- ggplot(pdat,aes(x = lfc1,y = lfc2)) +
  geom_point(color = "dodgerblue",shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,linetype = "dotted") +
  labs(x = "first posterior mean",y = "second posterior mean") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = z1,y = z2)) +
  geom_point(color = "dodgerblue",shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,linetype = "dotted") +
  labs(x = "first z-score estimate",y = "second z-score estimate") +
  xlim(-22,102) +
  ylim(-22,102) +
  theme_cowplot(font_size = 12)
plot_grid(p1,p2)

# Save the plots.
ggsave("../plots/mcmc_scatterplots_sims_k=6.png",
       plot_grid(p1,p2),height = 3,width = 6,dpi = 600)
