# A short script use to compile the susie GSEA results for the droplet
# data into a table.
library(Matrix)
library(susieR)
library(pathways)
source("../code/gsea.R")

# Load the gene sets.
data(gene_sets_mouse)

# Load the results of the gene set enrichment analysis.
load("../output/droplet/gsea-droplet-curated-only.RData")

# Extract the info for the gene sets that were included in the gene
# set enrichment analysis.
gene_set_info <- gene_sets_mouse$gene_set_info
rownames(gene_set_info) <- gene_set_info$id
gene_set_info <- gene_set_info[colnames(X),]

# Reorder the topics to correspond to the ordering used in the
# manuscript.
topics <- c("k2","k5","k7","k4","k1","k6","k3","k5+k7")
Y      <- Y[,topics]
gsea   <- gsea[topics]
colnames(Y) <- c(paste0("k",1:7),"k2+k3")
names(gsea) <- c(paste0("k",1:7),"k2+k3")

# Compile the enriched gene sets into a single table.
dat <- NULL
topics <- names(gsea)
for (i in topics) {
  x   <- compile_gsea_table(gsea[[i]],X,Y[,i],gene_set_info)
  x   <- cbind(data.frame(topic = i),x)
  dat <- rbind(dat,x)
}

# Write the data frame to a CSV file.
dat <-
  transform(dat,
            lbf = format(round(lbf,digits = 3),trim = TRUE,scientific = FALSE),
            pip = format(round(pip,digits = 3),trim = TRUE,scientific = FALSE),
            coef = format(round(coef,digits = 3),trim = TRUE,scientific=FALSE))
cat("Check for double quotes in description_brief field:",
    any(grepl("\"",dat$description_brief,fixed = TRUE)),"\n")
write.csv(dat,"gsea_droplet_curated_only.csv",quote = TRUE,row.names = FALSE)
