fit_cluster <- init_poisson_nmf_from_clustering(counts,samples$celltype)
fit_cluster <- poisson2multinom(fit_cluster)
loglik_topic <- loglik_multinom_topic_model(counts,fit)
loglik_cluster <- loglik_multinom_topic_model(counts,fit_cluster)

ggplot(data.frame(x = loglik_topic),aes(x)) +
  geom_histogram(color = "black",fill = "black",bins = 64) +
  theme_cowplot(font_size = 10) +
  labs(x = "loglik",y = "number of cells")

pdat <- data.frame(x = loglik_cluster,
                   y = loglik_topic,
                   maxp = apply(fit$L,1,max))
pdat <- subset(pdat,y > -5000)
ggplot(pdat,aes(x = x,y = y,fill = maxp)) +
  geom_point(shape = 21,color = "white") +
  scale_fill_gradient2(mid = "gold",high = "orangered",low = "deepskyblue",
                       midpoint = 0.5) +
  geom_abline(intercept = 0,slope = 1,linetype = "dotted",
              color = "black") +
  theme_cowplot(font_size = 12)
