# fit_clusters <- fit_multinom_model(samples$cluster,counts)

B <- le_lfc(fit_clusters$F,e = 1e-8)
D <- min_kl_poisson(fit_clusters$F)

# B cells
k   <- "B"
dat <- cbind(genes,data.frame(lfc = B[,k],kl = D[,k]))
print(subset(dat,lfc > 3 & kl > 0.001))

# CD14+ cells
k   <- "CD14+"
dat <- cbind(genes,data.frame(lfc = B[,k],kl = D[,k]))
print(subset(dat,lfc > 3 & kl > 0.001))

# CD34+ cells
k   <- "CD34+"
dat <- cbind(genes,data.frame(lfc = B[,k],kl = D[,k]))
print(subset(dat,lfc > 5 & kl > 3e-4))

# CD8+ T cells
k   <- "CD8+"
dat <- cbind(genes,data.frame(lfc = B[,k],kl = D[,k]))
print(subset(dat,lfc > 1.5 & kl > 1e-4))

# dendritic cells
k   <- "dendritic"
dat <- cbind(genes,data.frame(lfc = B[,k],kl = D[,k]))
print(subset(dat,lfc > 4 & kl > 1e-4))

# NK cells
k   <- "NK"
dat <- cbind(genes,data.frame(lfc = B[,k],kl = D[,k]))
print(subset(dat,lfc > 2 & kl > 0.001))
