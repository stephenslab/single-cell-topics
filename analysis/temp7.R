k <- "k5"
D <- fastTopics:::min_kl_poisson(fit$F + 1e-8)
D <- D[rows,]
pdat <- data.frame(gene     = genes$symbol,
                   lkl      = log10(D[,k] + 1e-10),
                   f0       = log10(de$f0),
                   postmean = de$postmean[,k],
                   z        = de$z[,k],
                   lfsr     = de$lfsr[,k],
                   stringsAsFactors = FALSE)
pdat <- subset(pdat,z > 0)
pdat$lfsr <- cut(pdat$lfsr,c(-1,0.001,0.01,0.05,Inf))
pdat <- transform(pdat,z = pmin(150,abs(z)))
pdat1 <- pdat
i <- which(pdat$z < 40)
pdat1[i,"gene"] <- ""
p1 <- ggplot(pdat1,aes(x = lkl,y = z,label = gene)) +
  geom_point() +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,na.rm = TRUE) +
  labs(x = "log10 K-L divergence",
       y = "|z-score|") +
  theme_cowplot()
pdat2 <- pdat
i <- which(!(pdat$postmean > 8 & pdat$z > 2))
pdat2[i,"gene"] <- ""
p2 <- ggplot(pdat2,aes(x = lkl,y = postmean,fill = lfsr,label = gene)) +
  geom_point(shape = 21,color = "white") +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,na.rm = TRUE) +
  scale_fill_manual(values = c("deepskyblue","gold","orange","coral")) +
  labs(x = "log10 K-L divergence",
       y = "posterior mean LFC") +
  theme_cowplot()
p3 <- ggplot(pdat,aes(x = lkl,fill = lfsr)) +
  geom_histogram(bins = 64) +
  scale_fill_manual(values = c("deepskyblue","gold","orange","coral")) +
  labs(x = "log10 K-L divergence",y = "genes") +
  theme_cowplot()
# p4 <- ggplot(pdat4,aes(x = -log10(lfsr + 1e-8),y = lkl)) +
#   geom_point(color = "white",fill = "royalblue",shape = 21) +
#   theme_cowplot()
pdat4 <- subset(pdat,lfsr > 0.05)
p4 <- ggplot(pdat4,aes(x = lkl,y = f0))+
  geom_point() +
  theme_cowplot()
