library(Matrix)
library(DESeq2)
library(MAST)
library(fastTopics)
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
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dashed") +
  labs(x = "first posterior mean",y = "second posterior mean") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = z1,y = z2)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dashed") +
  labs(x = "first z-score estimate",y = "second z-score estimate") +
  xlim(-45,42) + 
  ylim(-45,42) + 
  theme_cowplot(font_size = 12)

# Save the plots.
ggsave("../plots/mcmc_scatterplots_sims_k=2_alpha=0.01.png",
       plot_grid(p1,p2),height = 3,width = 6,dpi = 600)

# Compare DESeq2 and fastTopics output.
pdat1 <- data.frame(lfc.deseq      = combine_sim_res(res,function (x) x$deseq$log2FoldChange),
                    lfc.fasttopics = combine_sim_res(res,function (x) x$de1$postmean[,2]),
                    z.deseq        = combine_sim_res(res,
                                       function (x) with(x$deseq,log2FoldChange/lfcSE)),
                    z.fasttopics   = combine_sim_res(res,function (x) x$de1$z[,2]))
pdat2 <- data.frame(deseq2 = combine_sim_res(res,function (x) x$deseq$svalue),
                    shrink = combine_sim_res(res,function (x) x$de1$svalue[,2]),
                    de     = factor(combine_sim_res(res,
                               function (x) with(x$dat,abs(F[,1] - F[,2]) > 1e-8))))
p1 <- ggplot(pdat1,aes(x = lfc.deseq,y = lfc.fasttopics)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dashed") +
  labs(x = "DESeq2",y = "GoM DE",title = "LFC estimates") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat1,aes(x = z.deseq,y = z.fasttopics)) +
  geom_point(shape = 4,size = 0.75) +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dashed") +
  labs(x = "DESeq2",y = "GoM DE",title = "z-scores") +
  theme_cowplot(font_size = 12)
p3 <- ggplot(pdat2,aes(x = deseq2,color = de,fill = de)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  labs(x = "s-value",y = "genes",title = "DESeq2") +
  theme_cowplot(font_size = 12)
p4 <- ggplot(pdat2,aes(x = shrink,color = de,fill = de)) +
  geom_histogram(bins = 64,show.legend = FALSE) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  labs(x = "s-value",y = "genes",title = "fastTopics") +
  theme_cowplot(font_size = 12)
       
# Save the plots.
ggsave("../plots/deseq2_vs_fasttopics_sims.png",
       plot_grid(p1,p2,p3,p4,nrow = 2,ncol = 2),
       height = 3,width = 9,dpi = 600)
