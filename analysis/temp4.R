ns <- 10000
set.seed(1)
de1 <- de_analysis(fit,X,control = list(ns = ns,nc = 2))
set.seed(2)
de2 <- de_analysis(fit,X,control = list(ns = ns,nc = 2))
cor(de1$z[,2],de2$z[,2])
plot(de1$z[,2],de2$z[,2],pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
