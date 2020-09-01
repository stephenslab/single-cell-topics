# This function aligns the gene-set data (gene_info, gene_sets) with
# the gene-wise statistics (genes, diff_count_res) by their Ensembl
# ids to prepare these data for a gene-set enrichment analysis. It is
# assumed that gene_info and genes each have an "ensembl" column
# containing the Ensembl gene ids.
align_gene_data_by_ensembl <- function (gene_info, gene_sets, 
                                        genes, diff_count_res) {
  ids       <- intersect(gene_info$Ensembl,genes$ensembl)
  i         <- match(ids,gene_info$Ensembl)
  j         <- match(ids,genes$ensembl)
  gene_info <- gene_info[i,]
  gene_sets <- gene_sets[i,]
  genes     <- genes[j,]
  diff_count_res$colmeans <- diff_count_res$colmeans[j]
  diff_count_res$F0       <- diff_count_res$F0[j,]
  diff_count_res$F1       <- diff_count_res$F1[j,]
  diff_count_res$beta     <- diff_count_res$beta[j,]
  diff_count_res$se       <- diff_count_res$se[j,]
  diff_count_res$Z        <- diff_count_res$Z[j,]
  diff_count_res$pval     <- diff_count_res$pval[j,]
  return(list(gene_info = gene_info,gene_sets = gene_sets,
              genes = genes,diff_count_res = diff_count_res))
}
