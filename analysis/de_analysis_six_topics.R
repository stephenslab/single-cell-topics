# Having compared the overall distributions of the gene-wise rankings,
# we now compare performance of these rankings by plotting power
# vs. FDR. Again, we restrict this comparison only to positive LFCs.
v1 <- create_fdr_vs_power_curve(-pdat$lkl,nonzero_lfc,length.out = 400)
v2 <- create_fdr_vs_power_curve(pdat$noshrink,nonzero_lfc,length.out = 400)
v3 <- create_fdr_vs_power_curve(pdat$lfsr,nonzero_lfc,length.out = 400)
dat <- rbind(cbind(v1,method = "kl"),
             cbind(v2,method = "noshrink"),
             cbind(v3,method = "shrink"))
p <- ggplot(dat,aes(x = fdr,y = power,color = method)) +
  geom_line(size = 0.65,orientation = "y")  +
  scale_color_manual(values = c("royalblue","limegreen","orange")) +
  theme_cowplot()
print(p)

# Save the plots.
ggsave("../plots/fdr_vs_power_sims_k=6.eps",p,height = 2.5,width = 4)
