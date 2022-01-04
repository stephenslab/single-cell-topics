library(Matrix)
library(ggplot2)
library(cowplot)
source("../code/de_analysis_functions.R")

# Load the results of the simulations when simulated topic proportions
# are mostly 0 or 1.
load("../output/sims/sims-k=2-alpha=0.01.RData")
res["session.info"] <- NULL

# Most of the simulated cells should be attributed to a single topic.
combine_sim_res <- function (res, f)
  do.call("c",lapply(res,function (x) f(x)))
x <- combine_sim_res(res,function (x) apply(x$dat$L,1,max))
print(mean(x > 0.99))

# Assess accuracy of Monte Carlo estimates (for topic 2 only).
pdat <- data.frame(postmean1 = combine_sim_res(res,function (x) x$de1$postmean[,2]),
                   postmean2 = combine_sim_res(res,function (x) x$de2$postmean[,2]),
                   z1        = combine_sim_res(res,function (x) x$de1$z[,2]),
		   z2        = combine_sim_res(res,function (x) x$de2$z[,2]))
p1 <- ggplot(pdat,aes(x = postmean1,y = postmean2)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "darkorange",linetype = "dashed") +
  labs(x = "first posterior mean",y = "second posterior mean") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = z1,y = z2)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "darkorange",linetype = "dashed") +
  labs(x = "first z-score estimate",y = "second z-score estimate") +
  xlim(-45,42) + 
  ylim(-45,42) + 
  theme_cowplot(font_size = 12)

# Save the plots.
ggsave("../plots/mcmc_scatterplots_sims_k=2_alpha=0.01.png",
       plot_grid(p1,p2),height = 3,width = 6,dpi = 600)
