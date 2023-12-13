#procrustes and mantel 

data <- read.csv("asv.txt", head=TRUE,sep="\t",row.names = 1)
data <- data[which(rowSums(data) > 0),]
data <- t(data)
s.dist <- vegdist(data,method = "bray")

data <- read.csv("environment.txt", head=TRUE,sep="\t",row.names = 1)
r.dist <- vegdist(data)
#Mantel test
mantel(s.dist,r.dist)
#Spearman correlation analysis
mantel(s.dist,r.dist,method = "spearman")
#Dimensionality reduction
mds.s <- monoMDS(s.dist)
mds.r <- monoMDS(r.dist)
#Procrustes analysis
pro.s.r <- procrustes(mds.s,mds.r)
protest(mds.s,mds.r)

Y <- cbind(data.frame(pro.s.r$Yrot), data.frame(pro.s.r$X))
X <- data.frame(pro.s.r$rotation)
Y$ID <- rownames(Y)
#plot
p <- ggplot(Y) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#B2182B", size = 1) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#56B4E9", size = 1) +
  geom_point(aes(X1, X2), fill = "#B2182B", size = 4, shape = 21) +
  geom_point(aes(MDS1, MDS2), fill = "#56B4E9", size = 4, shape = 21) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Correlation between community and environment") + 
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = 'Procrustes analysis:\n    M2 = 0.8358, p-value = 0.035\nMantel test:\n    r = 0.1703, p-value = 0.04',
           x = -1.5, y = 1.2, size = 4,hjust = 0) +
  theme(plot.title = element_text(size=16,colour = "black",hjust = 0,face = "bold"))