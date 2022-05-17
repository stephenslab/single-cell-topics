# A short script to compile a CSV file containing the results of the
# grade-of-membership differential expression (GoM DE) analysis for
# the mixture of purified PBMC data.
library(fastTopics)

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
topics      <- c(3,2,5,4,1,6)
de$postmean <- de$postmean[,topics]
de$z        <- de$z[,topics]
de$lfsr     <- de$lfsr[,topics]
    
# Compile the DE results for all topics into a single table.
dat <- NULL
k <- ncol(de$z)
for (i in 1:k) {
  x <- data.frame(topic   = i,
                  gene    = genes$symbol,
                  ensembl = genes$ensembl,
                  lfc     = de$postmean[,i],
                  z       = de$z[,i],
                  lfsr    = de$lfsr[,i],
                  stringsAsFactors = FALSE)
  dat <- rbind(dat,x)
}

# Filter out genes with lfsr >= 0.01.
dat <- subset(dat,lfsr < 0.01)

# Write the data frame to a CSV file. Here we have
dat <- transform(dat,
                 topic = factor(topic),
                 lfc   = round(lfc,digits = 3),
                 z     = round(z,digits = 3),
                 lfsr  = format(lfsr,digits = 3,trim = TRUE))
write.csv(dat,"de_purified_pbmc.csv",quote = FALSE,row.names = FALSE)
