####Figure2 the beta diversity
df<-read.csv("asv.csv")
#pcoa
df.pca<-prcomp(df[,1:7],scale. = T)
pca.results<-data.frame(df.pca$x)[,1:2]
pca.results$target<-paste0('cultivar',df$target)
#Computing center point
centroid <- aggregate(cbind(PC1,PC2) ~ target,
                      data = pca.results,
                      FUN = mean)
#Combined with the results of pcoa
pca.results1<-dplyr::left_join(pca.results, centroid, by = "target", 
                               suffix = c("",".cen"))
cols=c("#B3808B","#537874")
ggplot(pca.results1,aes(x=PC1,y=PC2))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=target))+
  geom_segment(aes(xend=PC1.cen,yend=PC2.cen,color=target),
               show.legend = F)+
  ggtheme+
  scale_colour_manual(values = cols)+
  geom_label(data = centroid, 
             aes(label = target, fill = target), size = 5, 
             show.legend = FALSE,
             color="white")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "top")

#Bray-Curtis distances 
dis <- vegdist(varespec)
#groups
groups <- factor(c(rep(1,99)), labels = c("MPs","Soil"))
#Calculate multivariate dispersions
mod <- betadisper(dis, groups)
#ANOVA
anova(mod)
#Permutation test
permutest(mod, pairwise = TRUE, permutations = 99)
#Tukey's Honest Significant Differences
mod.HSD <- TukeyHSD(mod)
boxplot(mod)

