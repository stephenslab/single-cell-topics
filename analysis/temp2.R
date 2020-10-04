rows1 <- which(with(samples_droplet,
                    cluster == "T+N" & tissue == "Neuroendocrine"))
rows2 <- which(with(samples_droplet,
                    cluster == "T+N" & tissue == "Tuft"))
summary(poisson2multinom(fit_droplet)$L[rows1,c(2,6)])
summary(poisson2multinom(fit_droplet)$L[rows2,c(2,6)])
