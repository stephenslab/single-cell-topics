samples <- readRDS("../output/pbmc-68k/clustering-pbmc-68k.rds")
fit_merged <- merge_topics(fit,c("k3","k5","k10"))
fit_merged <- merge_topics(fit_merged,c("k2","k7"))
fit_merged <- merge_topics(fit_merged,c("k9","k11"))
fit_merged <- merge_topics(fit_merged,c("k6","k8"))
rows <- which(samples$cluster == "CD8+" |
              samples$cluster == "TA" |
              samples$cluster == "TB" |
              samples$cluster == "NK")
fit_t <- select_loadings(fit,loadings = rows)
# out1 <- diff_count_analysis(fit_merged,counts)
# out2 <- diff_count_analysis(fit,counts)
# out3 <- diff_count_clusters(samples$cluster,counts)
diff_count_res <- diff_count_analysis(fit_t,counts[rows,])
# volcano_plot(out3,k = "CD8+",labels = genes$symbol,
#              label_above_quantile = 0.998)

source("../unused/gsea.R")
genesetfile  <- "../data/gene_sets_human.RData"
load("~/git/pathways/output/gene_sets_human.RData")
rownames(gene_sets) <- gene_info$Ensembl

# Prepare the gene-set data and gene-wise statistics for the gene-set
# enrichment analysis. First, align the gene-set data with the
# gene-wise statistics.
out            <- align_gene_data(gene_sets,diff_count_res)
gene_sets      <- out$gene_sets
diff_count_res <- out$diff_count_res
ids            <- rownames(gene_sets)
gene_info      <- gene_info[match(ids,gene_info$Ensembl),]
genes          <- genes[match(ids,genes$ensembl),]
rm(out,ids)

# Next, remove gene sets with fewer than 4 genes, and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(gene_sets) >= 4 & colSums(gene_sets) <= 400)
gene_set_info <- gene_set_info[i,]
gene_sets     <- gene_sets[,i]
rm(i)

# For each topic, perform a gene-set enrichment analysis using fgsea.
gsea_res <- perform_gsea_all_topics(gene_sets,diff_count_res,nproc = 4)
gsea_res$ES[is.na(gsea_res$ES)] <- 0
gsea_res$NES[is.na(gsea_res$NES)] <- 0
gsea_res$pval[is.na(gsea_res$pval)] <- 1

# Create a plotly scatterplot to explore the gene-set enrichment
# results.
p1 <- create_gsea_plotly(gene_set_info,gsea_res,5)
saveWidget(p1,"gsea.html",selfcontained = TRUE)
