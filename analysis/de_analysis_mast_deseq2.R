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
