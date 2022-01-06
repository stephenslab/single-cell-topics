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

# Next we consider the K-L divergence measure used in Dey, Hsiao &
# Stephens (2017) to rank genes, and compare its ranking to a ranking
# based on *p*-values (without using adaptive shrinkage) or *s*-values
# (after applying adaptive shrinkage). Since the K-L divergence is not
# a signed ranking, we restrict this comparison only to LFCs that are
# estimated to be positive.
get_nonneg_lfcs <- function (de)
  which(de$est > -1e-8)
nonzero_lfc <- combine_sim_res(res,
                 function (x) {
		   n <- nrow(x$dat$F)
		   k <- ncol(x$dat$F)
		   out <- matrix(FALSE,n,k)
                   for (i in 1:n) {
                     y <- x$dat$F[i,]
                     if (max(y) - min(y) > 1e-8) {
                       j <- which.max((y - mean(y))^2)
                       out[i,j] <- TRUE
                     }
		   }
		   return(out[get_nonneg_lfcs(x$de1)])
		 })
noshrink <-
  combine_sim_res(res,function (x) x$de0$lpval[get_nonneg_lfcs(x$de1)])
shrink <-
  combine_sim_res(res,function(x)x$de1$svalue[get_nonneg_lfcs(x$de1)])
lfsr <- combine_sim_res(res,function(x)x$de1$lfsr[get_nonneg_lfcs(x$de1)])
lkl <- combine_sim_res(res,function (x)
         fastTopics:::min_kl_poisson(x$fit$F)[get_nonneg_lfcs(x$de1)])
pdat <- data.frame(nonzero_lfc = factor(nonzero_lfc),
                   noshrink    = 10^(-noshrink),
                   shrink      = shrink,
                   lfsr        = lfsr,
                   lkl         = log10(lkl + 1e-8))
p1 <- ggplot(pdat,aes(x = lkl,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  coord_cartesian(ylim = c(0,2.5e4)) +
  labs(x = "log10 K-L divergence",y = "genes",title = "K-L divergence") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = noshrink,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  coord_cartesian(ylim = c(0,2.5e4)) +
  labs(x = "p-value",y = "genes",title = "without shrinkage") +
  theme_cowplot(font_size = 12)
p3 <- ggplot(pdat,aes(x = shrink,color = nonzero_lfc,fill = nonzero_lfc)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  coord_cartesian(ylim = c(0,2.5e4)) +
  labs(x = "s-value",y = "genes",title = "with shrinkage") +
  theme_cowplot(font_size = 12)
plot_grid(p1,p2,p3,nrow = 1,ncol = 3)

# Save the plots.
ggsave("../plots/rankings_sims_k=6.eps",
       plot_grid(p1,p2,p3,nrow = 1,ncol = 3),height = 2.25,width = 7)

# Having compared the overall distributions of the gene-wise rankings,
# we now compare performance of these rankings by plotting power
# vs. FDR. Again, we restrict this comparison only to positive LFCs.
v1  <- create_fdr_vs_power_curve(-pdat$lkl,nonzero_lfc,length.out = 400)
v2  <- create_fdr_vs_power_curve(pdat$lfsr,nonzero_lfc,length.out = 400)
dat <- rbind(cbind(v1,method = "kl"),
             cbind(v2,method = "fastTopics"))
p <- ggplot(dat,aes(x = fdr,y = power,color = method)) +
  geom_line(size = 0.65,orientation = "y")  +
  scale_color_manual(values = c("royalblue","orange")) +
  theme_cowplot()
print(p)

# Save the plots.
ggsave("../plots/fdr_vs_power_sims_k=6.eps",p,height = 2.5,width = 4)
