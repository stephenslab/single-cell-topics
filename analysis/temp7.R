de             <- de1[c("f0","est","postmean","z","lfsr")]
class(de)      <- c("topic_model_de_analysis","list")
i              <- which(abs(de2$z) < abs(de1$z))
de$postmean[i] <- de2$postmean[i]
de$z[i]        <- de2$z[i]
de$lfsr[i]     <- de$lfsr[i]

rows        <- match(rownames(deseq$z),rownames(de$z))
de$f0       <- de$f0[rows]
de$est      <- de$est[rows,]
de$postmean <- de$postmean[rows,]
de$z        <- de$z[rows,]
de$lfsr     <- de$lfsr[rows,]

D <- fastTopics:::min_kl_poisson(fit$F + 1e-8)
D <- D[rows,]

k <- "k3"

pdat <- data.frame(lkl      = log10(D[,k] + 1e-10),
                   est      = de$est[,k],
                   postmean = de$postmean[,k],
                   lfsr     = de$lfsr[,k],
                   stringsAsFactors = FALSE)
pdat <- subset(pdat,est > 0)
pdat <- transform(pdat,lfsr = cut(lfsr,c(-1,0.001,0.01,0.05,Inf)))

p1 <- ggplot(pdat,aes(x = lkl,y = est,fill = postmean)) +
  geom_point(shape = 21,color = "white") +
  scale_fill_gradient2(low = "limegreen",mid = "gold",high = "orangered",
                       midpoint = 4) +
  labs(x = "log10 K-L divergence",y = "MLE") +
  theme_cowplot(font_size = 12)

p2 <- ggplot(pdat,aes(x = lkl,y = est,fill = lfsr)) +
  geom_point(shape = 21,color = "white") +
  scale_fill_manual(values = c("deepskyblue","gold","orange","coral")) +
  labs(x = "log10 K-L divergence",
       y = "MLE") +
  theme_cowplot(font_size = 12)
plot_grid(p1,p2)
