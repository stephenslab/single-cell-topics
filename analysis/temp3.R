j <- which(is.element(genes$symbol,rpgenes))
y <- rowSums(counts[,j])/rowSums(counts)
ggplot(data.frame(x = samples$celltype,y = y),aes(x = x,y = y)) +
  geom_boxplot(width = 0.5) +
  labs(x = "",y = "proportion") +
  theme_cowplot(font_size = 10) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
    
