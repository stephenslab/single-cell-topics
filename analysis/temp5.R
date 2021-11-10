library(Seurat)
library(MAST)
seurat <- CreateSeuratObject(counts = t(X))
Idents(seurat) <- cluster
mast <- FindMarkers(seurat,ident.1 = "2",ident.2 = NULL,test.use = "MAST",
                    logfc.threshold = 0,min.pct = 0)

# Compare the LFC estimates.
genes <- rownames(mast)
pdat  <- data.frame(MAST       = mast$avg_logFC,
                    fastTopics = de$postmean[genes,2])
ggplot(pdat,aes(x = MAST,y = fastTopics)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "skyblue",linetype = "dotted") +
  ggtitle("LFC estimates") +
  theme_cowplot()

# Compare the p-values.
pdat <- data.frame(MAST       = pmax(-10,log10(mast$p_val)),
                   fastTopics = pmax(-10,-de$lpval[genes,2]))
ggplot(pdat,aes(x = MAST,y = fastTopics)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "skyblue",linetype = "dotted") +
  ggtitle("p-values") +
  theme_cowplot()

# Plot power vs. FDR.
mast_pval                 <- rep(1,10000)
names(mast_pval)          <- colnames(X)
mast_pval[rownames(mast)] <- mast$p_val
pdat1 <- create_fdr_vs_power_curve(deseq$padj,nonzero_lfc,length.out = 200)
pdat2 <- create_fdr_vs_power_curve(de$lfsr[,2],nonzero_lfc,length.out = 200)
pdat3 <- create_fdr_vs_power_curve(mast_pval,nonzero_lfc,length.out = 200)
pdat  <- rbind(cbind(pdat1,method = "DESeq2"),
               cbind(pdat2,method = "fastTopics"),
               cbind(pdat3,method = "MAST"))
ggplot(pdat,aes(x = fdr,y = power,color = method)) +
  geom_line(size = 0.65,orientation = "y")  +
  scale_color_manual(values = c("darkblue","tomato","dodgerblue")) +
  theme_cowplot()

# Plot ROC curve.
pdat1 <- create_roc_curve(deseq$padj,nonzero_lfc,length.out = 100)
pdat2 <- create_roc_curve(de$lfsr[,2],nonzero_lfc,length.out = 100)
pdat3 <- create_roc_curve(mast_pval,nonzero_lfc,length.out = 100)
pdat  <- rbind(cbind(pdat1,method = "DESeq2"),
               cbind(pdat2,method = "fastTopics"),
               cbind(pdat3,method = "MAST"))
ggplot(pdat,aes(x = fpr,y = tpr,color = method)) +
  geom_line(size = 0.65)  +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  scale_color_manual(values = c("darkblue","tomato","dodgerblue")) +
  theme_cowplot()
