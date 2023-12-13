###figure 3 the distance-decay relationship
asv <- read.csv("genus_relative.txt", sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
colnames(asv) <- asv[1,] 
for(i in 1: ncol(asv)){
  
  asv[,i]<-as.numeric(asv[,i])
  
}
#tree
roottree<- read_tree('tree.nwk')
#build the phyloseq object
physeq <- phyloseq(otu_table(asv, taxa_are_rows = TRUE), phy_tree(roottree))
#Calculate the weighted Unifrac distance
wunifrac<- distance(physeq, method = 'wunifrac')
#Calculate the Bray-Curtis Unifrac distance
bray <- distance(physeq, method = 'bray')
#time distance
distance <- read.csv("distance.txt", header = TRUE, sep = '\t', row.names = 1)
dist_distance <- vegdist(distance,method = "euclidean")
dist_distance <- log(dist_distance)
distance1 <- as.data.frame(as.vector(dist_distance))
colnames(data) <- c("d_distance","wunifrac","bray")
#rbind
data <- data.frame(as.vector(distance1),as.vector(1-wunifrac))
summary(lm(data$bray~data$d_distance))
#plot
data.long <- pivot_longer(data, cols = -d_time, names_to = "type", values_to = "dissimilarity")
data.long$type <- fct_inorder(data.long$type)

ggplot(data.long, aes(x = d_distance,y = dissimilarity)) +
  geom_point() +
  geom_smooth(method = "lm",alpha = 0.2) +
  labs(x = "Distance (1000 Km)", y = "Community similarity (%) ") + theme_bw()
