# TO DO: Explain here what this function does, and how to use it.
create_abundance_plot <- function (fit, probs = c(0.1,0.25,0.5),
                                   clrs=c("tomato","dodgerblue","darkblue")) {
  m <- length(probs)
  k <- ncol(fit$L)
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  L      <- fit$L
  out    <- apply(L,2,function (x) rev(cumsum(rev(table(cut(x,c(probs,1)))))))
  topics <- colnames(L)
  topics <- topics[order(out[1,],decreasing = TRUE)]
  dat    <- data.frame(topic = factor(rep(colnames(L),each = m),topics),
                       prob  = factor(rep(probs,times = k)),
                       n     = as.vector(out))
  return(ggplot(dat,aes_string(x = "topic",y = "n",fill = "prob")) +
         geom_col(color = "white",position = "dodge",width = 0.8) +
         scale_fill_manual(values = clrs) +
         labs(x = "topic",y = "samples") +
         theme_cowplot(font_size = 10))
}
