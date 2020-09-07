p1 <- volcano_plot(diff_count_purified,k = 4,labels = genes_purified$symbol)

# Examine in greater detail relationship between topic 4 proportions
# and expression level of NKG7 in purified data.
load("../data/pbmc_purified.RData")
counts <- counts[,genes_purified$ensembl]
rm(samples,genes)
i    <- which(genes_purified$symbol == "NKG7")
fit2 <- poisson2multinom(fit_purified)
table(cut(fit2$L[,"k4"],seq(0,1,0.1)))
table(cut(fit2$L[samples_purified$cluster == "A1","k4"],seq(0,1,0.1)))
table(cut(fit2$L[samples_purified$cluster != "A1","k4"],seq(0,1,0.1)))
pdat <- data.frame(x = cut(fit2$L[,"k4"],c(-1,seq(0.1,1,0.1))),
                   y = counts[,i])
p3 <- ggplot(pdat,aes(x = x,y = y)) +
  geom_boxplot(width = 0.25,size = 0.4,outlier.shape = NA) +
  theme_cowplot(font_size = 10)

# Examine in greater detail relationship between topic 3 proportions
# and expression level of NKG7 in 68k data.
load("../data/pbmc_68k.RData")
counts <- counts[,genes_68k$ensembl]
rm(samples,genes)
i    <- which(genes_68k$symbol == "NKG7")
fit2 <- poisson2multinom(fit_68k)
table(cut(fit2$L[,"k3"],seq(0,1,0.1)))
pdat <- data.frame(x = cut(fit2$L[,"k3"],
                           c(-1,seq(0.1,1,0.1))),
                   y = counts[,i])
p4 <- ggplot(pdat,aes(x = x,y = y)) +
  geom_boxplot(width = 0.25,size = 0.4,outlier.shape = NA) +
  theme_cowplot(font_size = 10)
