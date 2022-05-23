# A short script to compile a CSV file containing the results of the
# grade-of-membership differential expression (GoM DE) analysis for
# the mixture of purified PBMC data.
library(fastTopics)
source("../code/de.R")

# Load the results of the topic-model-based DE analysis.
load("../output/pbmc-purified/de-pbmc-purified-seed=1.RData")
de1 <- de
load("../output/pbmc-purified/de-pbmc-purified-seed=2.RData")
de2 <- de
rm(de)

# When the two z-scores disagree, use the one that is nearer to zero.
de             <- de1[c("postmean","z","lfsr")]
i              <- which(abs(de2$z) < abs(de1$z))
de$postmean[i] <- de2$postmean[i]
de$z[i]        <- de2$z[i]
de$lfsr[i]     <- de$lfsr[i]

# Reorder the topics to correspond to the ordering used in the
# manuscript.
topics      <- c("k3","k2","k5","k4","k1","k6")
de$postmean <- de$postmean[,topics]
de$z        <- de$z[,topics]
de$lfsr     <- de$lfsr[,topics]
colnames(de$postmean) <- paste0("k",1:6)
colnames(de$z)        <- paste0("k",1:6)
colnames(de$lfsr)     <- paste0("k",1:6)

# Compile the DE results for all topics into a single table.
dat <- compile_de_table(de)

# Add the gene symbols.
ids <- dat$gene
rownames(genes) <- genes$ensembl
dat <- cbind(dat["topic"],genes[ids,],dat[c("lfc","z","lfsr")])

# Filter out genes with lfsr >= 0.01.
dat <- subset(dat,lfsr < 0.01)

# Reorder the genes by topic, then by LFC.
rows <- with(dat,order(topic,-lfc))
dat  <- dat[rows,]

# Write the data frame to a CSV file.
dat <-
  transform(dat,
            lfc = format(round(lfc,digits = 3),trim = TRUE,scientific = FALSE),
            z = format(round(z,digits = 3),trim = TRUE,scientific = FALSE),
            lfsr = format(lfsr,digits = 3,trim = TRUE,scientific = TRUE))
write.csv(dat,"de_purified_pbmc.csv",quote = FALSE,row.names = FALSE)
