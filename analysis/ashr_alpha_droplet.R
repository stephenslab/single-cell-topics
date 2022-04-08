library(fastTopics)
library(ashr)
set.seed(1)

# Load the results of the DE analysis.
load("../output/droplet/de-droplet-noshrink.RData")
b  <- c(de_merged$postmean)
se <- c(with(de_merged,postmean/z))

# Run adaptive shrinkage for different settings of alpha.
a <- seq(0,1,0.05)
n <- length(a)
fits <- vector("list",n)
for (i in 1:n) {
  cat(i,"")
  fits[[i]] <- ash(b,se,alpha = a[i])
}
cat("\n")
