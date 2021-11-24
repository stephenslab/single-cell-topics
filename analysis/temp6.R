# T cells
i <- "T cell"
k <- "k1"
pdat <- data.frame(gene                = genes$symbol,
                   postmean.deseq      = deseq$postmean[,i],
                   postmean.fasttopics = de$postmean[,k],
		   z.deseq             = deseq$z[,i],
		   z.fasttopics        = de$z[,k],
		   lfsr = cut(de$lfsr[,k],c(-1,0.001,0.01,0.05,Inf)),
		   stringsAsFactors = FALSE)
j <- which(with(pdat,
                !(postmean.fasttopics > 6 |
                  (postmean.fasttopics > 4 & postmean.deseq < 3))))
pdat[j,"gene"] <- ""
pdat <- subset(pdat,abs(z.deseq) > 2 | abs(z.fasttopics) > 2)
ggplot(pdat,aes(x = postmean.deseq,y = postmean.fasttopics,
                fill = lfsr,label = gene)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,na.rm = TRUE) +
  scale_fill_manual(values = c("deepskyblue","gold","orange","coral"),
                    na.value = "gainsboro") +
  labs(x = "DESeq2",y = "fastTopics") +
  theme_cowplot()

# NK cells
i <- "CD56+ NK"
k <- "k4"
pdat <- data.frame(gene                = genes$symbol,
                   postmean.deseq      = deseq$postmean[,i],
                   postmean.fasttopics = de$postmean[,k],
		   z.deseq             = deseq$z[,i],
		   z.fasttopics        = de$z[,k],
		   lfsr = cut(de$lfsr[,k],c(-1,0.001,0.01,0.05,Inf)),
		   stringsAsFactors = FALSE)
j <- which(with(pdat,
                !(postmean.fasttopics > 8 |
                  (postmean.fasttopics > 4 & postmean.deseq < 2))))
pdat[j,"gene"] <- ""
pdat <- subset(pdat,abs(z.deseq) > 2 | abs(z.fasttopics) > 2)
ggplot(pdat,aes(x = postmean.deseq,y = postmean.fasttopics,
                fill = lfsr,label = gene)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,na.rm = TRUE) +
  scale_fill_manual(values = c("deepskyblue","gold","orange","coral"),
                    na.value = "gainsboro") +
  labs(x = "DESeq2",y = "fastTopics") +
  theme_cowplot()

# CD14+ cells
i <- "CD14+ Monocyte"
k <- "k2"
pdat <- data.frame(gene                = genes$symbol,
                   postmean.deseq      = deseq$postmean[,i],
                   postmean.fasttopics = de$postmean[,k],
		   z.deseq             = deseq$z[,i],
		   z.fasttopics        = de$z[,k],
		   lfsr = cut(de$lfsr[,k],c(-1,0.001,0.01,0.05,Inf)),
		   stringsAsFactors = FALSE)
j <- which(with(pdat,
                !(postmean.fasttopics > 10 |
                  (postmean.fasttopics > 5 & postmean.deseq < 2))))
pdat[j,"gene"] <- ""
pdat <- subset(pdat,abs(z.deseq) > 2 | abs(z.fasttopics) > 2)
ggplot(pdat,aes(x = postmean.deseq,y = postmean.fasttopics,
                fill = lfsr,label = gene)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,na.rm = TRUE) +
  scale_fill_manual(values = c("deepskyblue","gold","orange","coral"),
                    na.value = "gainsboro") +
  labs(x = "DESeq2",y = "fastTopics") +
  theme_cowplot()

# CD34+ cells
i <- "CD34+"
k <- "k5"
pdat <- data.frame(gene                = genes$symbol,
                   postmean.deseq      = deseq$postmean[,i],
                   postmean.fasttopics = de$postmean[,k],
		   z.deseq             = deseq$z[,i],
		   z.fasttopics        = de$z[,k],
		   lfsr = cut(de$lfsr[,k],c(-1,0.001,0.01,0.05,Inf)),
		   stringsAsFactors = FALSE)
j <- which(with(pdat,postmean.fasttopics < 9))
pdat[j,"gene"] <- ""
pdat <- subset(pdat,abs(z.deseq) > 2 | abs(z.fasttopics) > 2)
ggplot(pdat,aes(x = postmean.deseq,y = postmean.fasttopics,
                fill = lfsr,label = gene)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,na.rm = TRUE) +
  scale_fill_manual(values = c("deepskyblue","gold","orange","coral"),
                    na.value = "gainsboro") +
  labs(x = "DESeq2",y = "fastTopics") +
  theme_cowplot()

# Ribosomal protein genes
k <- "k6"
pdat <- data.frame(gene     = genes$symbol,
                   postmean = de$postmean[,k],
                   z        = pmin(200,abs(de$z[,k])),
                   lfsr     = cut(de$lfsr[,k],c(-1,0.001,0.01,0.05,Inf)),
                   stringsAsFactors = FALSE)
j <- which(with(pdat,!(postmean > 4 | z > 100)))
pdat[j,"gene"] <- ""
ggplot(pdat,aes(x = postmean,y = z,fill = lfsr,label = gene)) +
  geom_point(shape = 21,color = "white",size = 1.5) +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,na.rm = TRUE) +
  scale_y_continuous(trans = "sqrt",breaks = c(1,2,5,10,20,50,100)) +
  scale_fill_manual(values = c("deepskyblue","gold","orange","tomato"),
                    na.value = "gainsboro") +
  labs(x = "posterior mean LFC",y = "|z-score|") +
  theme_cowplot()
