power1 <- create_fdr_vs_power_curve(-D[i],nonzero_lfc[i],length.out = 400)
power2 <- create_fdr_vs_power_curve(de$lfsr[i],nonzero_lfc[i],length.out=400)
roc1 <- create_roc_curve(-log(D[i]),nonzero_lfc[i],length.out = 400)
roc2 <- create_roc_curve(de$lfsr[i],nonzero_lfc[i],length.out = 400)
pdat1  <- cbind(power1,roc1[c("tpr","fpr")])
pdat2  <- cbind(power2,roc2[c("tpr","fpr")])
p1 <- ggplot(pdat1,aes(x = fdr,y = power,color = tpr)) +
  geom_point(size = 2)  +
  scale_color_gradient2(low = "deepskyblue",mid = "gold",
                        high = "orangered",midpoint = 0.5) +
  ggtitle("kl") +
  theme_cowplot()
p2 <- ggplot(pdat1,aes(x = fdr,y = power,color = fpr)) +
  geom_point(size = 2)  +
  scale_color_gradient2(low = "deepskyblue",mid = "gold",
                        high = "orangered",midpoint = 0.5) +
  ggtitle("kl") +
  theme_cowplot()
p3 <- ggplot(pdat2,aes(x = fdr,y = power,color = tpr)) +
  geom_point(size = 2)  +
  scale_color_gradient2(low = "deepskyblue",mid = "gold",
                        high = "orangered",midpoint = 0.5) +
  ggtitle("fastTopics") +
  theme_cowplot()
p4 <- ggplot(pdat2,aes(x = fdr,y = power,color = fpr)) +
  geom_point(size = 2)  +
  scale_color_gradient2(low = "deepskyblue",mid = "gold",
                        high = "orangered",midpoint = 0.5) +
  ggtitle("fastTopics") +
  theme_cowplot()
print(plot_grid(p1,p2,p3,p4))
