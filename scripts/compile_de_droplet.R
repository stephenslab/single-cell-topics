# A short script to compile a CSV file containing the results of the
# grade-of-membership differential expression (GoM DE) analysis for
# the droplet data.
library(fastTopics)
source("../code/de.R")

# Load the results of the topic-model-based DE analysis.
load("../output/droplet/de-droplet.RData")
out <- readRDS("../output/droplet/refit-droplet.rds")
de_rare <- out$de
de_rare_merged <- out$de_merged

# Compile the DE results for all topics into a single table,
# reorderinng the topics to correspond to the ordering used in the
# manuscript.
dat <- rbind(compile_de_table(de_merged),
             compile_de_table(de,c("k5","k7")))
dat_rare <- rbind(compile_de_table(de_rare_merged),
                  compile_de_table(de_rare,c("k3","k4")))
levels(dat$topic) <- c("k5","k1","k7","k4","k2+k3","k6","k2","k3")
levels(dat_rare$topic) <- c("k11","k10","k8+k9","k12","k9","k8")
dat <- rbind(dat,dat_rare)
topics <- c(paste0("k",1:12),c("k2+k3","k8+k9"))
dat$topic <- factor(dat$topic,topics)

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
write.csv(dat,"de_droplet.csv",quote = FALSE,row.names = FALSE)
