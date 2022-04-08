library(fastTopics)
library(ashr)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load the results of the DE analysis.
load("../output/droplet/de-droplet-noshrink.RData")
b  <- c(de_merged$postmean)
se <- c(with(de_merged,postmean/z))
i  <- which(!is.na(se))
b  <- b[i]
se <- se[i]

# Run adaptive shrinkage for different settings of alpha.
a <- seq(0,1,0.05)
n <- length(a)
fits <- vector("list",n)
for (i in 1:n) {
  cat(i,"")
  fits[[i]] <- ash(b,se,alpha = a[i])
}
cat("\n")

# Plot likelihood vs. alpha.
loglik <- sapply(fits,"[[","loglik")
pdat1 <- data.frame(alpha = a,loglik = loglik)
p1 <- ggplot(pdat1,aes(x = alpha,y = loglik)) +
  geom_point() +
  geom_line() +
  ggtitle("loglik vs. alpha") +
  theme_cowplot(font_size = 12)

# Plot posterior estimates for best fit vs. worst fit.
pdat2 <- data.frame(b1 = fits[[1]]$result$PosteriorMean,
                    b2 = fits[[n]]$result$PosteriorMean,
                    se = se)
p2 <- ggplot(pdat2,aes(x = b1,y = b2,color = se)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  scale_color_gradient2(low = "lightskyblue",mid = "gold",high = "orangered",
                        midpoint = 3) +
  labs(x = "alpha = 0",y = "alpha = 1",color = "s.e.",
       title = "posterior estimates") +
  theme_cowplot(font_size = 12)

# Compare the distributions of the estimated coefs.
pdat3 <-
  data.frame(q1 = quantile(fits[[1]]$result$PosteriorMean,seq(0.05,0.95,0.01)),
             q2 = quantile(fits[[n]]$result$PosteriorMean,seq(0.05,0.95,0.01)))
x <- with(pdat3,c(q1,q2))
p3 <- ggplot(pdat3,aes(x = q1,y = q2)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "skyblue",linetype = "dashed") +
  coord_cartesian(xlim = range(x),ylim = range(x)) +
  labs(x = "alpha = 0",y = "alpha = 1",title = "quantiles") +
  theme_cowplot(font_size = 12)
plot_grid(p1,p2,p3,nrow = 1,ncol = 3)
