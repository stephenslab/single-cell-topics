# Droplet.
loadings_plot(poisson2multinom(fit),samples$tissue)
pca_plot(poisson2multinom(fit),pcs = 1:2)

# Pulse-seq.
sum(loglik_poisson_nmf(counts,fit))
sum(loglik_poisson_nmf(counts,fits[["fit-pulseseq-scd-ex-k=13"]]))
loadings_plot(poisson2multinom(fit),samples$tissue)
