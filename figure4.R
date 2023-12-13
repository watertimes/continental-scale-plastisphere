##figure 4 the null model for the microbial assembly 
#data Loading
otu <- read.csv('asv.csv', row.names = 1)
otu <- data.frame(t(otu))
#tree
library(ape)
tree <- read.tree('tree.nwk')
##Calculate MNTD and NTI using the picante package
library(picante)
#make sure the names on the phylogeny are ordered the same as the names in otu table
match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);
#calculate empirical betaMNTD
beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'spe.MNTD_weighted.csv',quote=F);
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE
#calculate randomized betaMNTD
beta.reps = 999; # number of randomizations
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);
for (rep in 1:beta.reps) {  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));  
  print(c(date(),rep));  
}
weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);
for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");    
  };
};
rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"spe.weighted_bNTI.csv",quote=F);
pdf("spe.weighted_bNTI_Histogram.pdf")
hist(weighted.bNTI)
dev.off()

library(permute);library(gee);library(vegan);
library(ape);library(picante);library(ecodist);
#########
######## calculate Bray-Curtis########
#########
library(ecodist)
bray.out = as.matrix(distance(otu,method = 'bray-curtis'))
write.csv(bray.out,"bray_weighted.csv");
####################
########## for Raup-Crick ##########
####################
metric = 'RC'
abund.for.names = 'weighted'
colnames(otu) = gsub(pattern = "X",replacement = "",x = colnames(otu))
print(dim(otu))
no.of.samples = nrow(otu)
metric = 'RC'
abund = T
abund.for.names = 'weighted'
source("Raup_Crick_Abundance_One_Comparison.r")
for (i in 1:(no.of.samples - 1)){
  for (j in (i + 1):(no.of.samples)){
    raup.crick.out = raup_crick_abundance_one_comparison(null.one.use = i,null.two.use = j,otu, plot_names_in_col1=FALSE, classic_metric=FALSE, split_ties=TRUE, reps=999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE);
    print(raup.crick.out)
    write.csv(raup.crick.out,paste('RC_',abund.for.names,'_Diss_comm1_',i,'_comm2_',j,'.csv',sep=""),row.names=T,quote=F);
  }
}
## find RC reps that still need to be done
all.files = list.files(pattern = 'RC_')
all.RC = all.files[grep("RC_weighted_Diss_comm1_",all.files)]
length(all.RC)
bray = read.csv('bray_weighted.csv',row.names=1);
raup.crick.out = bray
raup.crick.out[,] = -999
for (i in 1:length(all.RC)) {
  RC.temp = read.csv(paste(all.RC[i],sep=""),row.names=1)
  if (nrow(RC.temp) != 1 ) {print(c(i,"ERROR"))} 
  if (ncol(RC.temp) != 1 ) {print(c(i,"ERROR"))} 
  raup.crick.out[which(rownames(raup.crick.out) == rownames(RC.temp)),which(colnames(raup.crick.out) == colnames(RC.temp))] = RC.temp[1,1]
  print(RC.temp)
}
raup.crick.out = as.data.frame(as.matrix(as.dist(raup.crick.out)))
head(raup.crick.out)
raup.crick.out[1:5,1:5]
write.csv(raup.crick.out,'RC_weighted.csv',quote=F);
pdf(paste('RC.',abund.for.names,'_Histogram.pdf',sep=""))
raup.crick.out.for.hist = as.dist(raup.crick.out)
raup.crick.out.for.hist = raup.crick.out.for.hist[raup.crick.out.for.hist != -999]
hist(raup.crick.out.for.hist,breaks=40)
abline(v=c(-0.95,0.95),col=2,lwd=2,lty=2)
dev.off()

bNTI<-read.table("spe.weighted_bNTI.csv",header=T,sep=",",row.names=1)
rc<-read.table("RC_weighted.csv",header=T,sep=",",row.names=1)
#influences of different processes
## match names in bNTI and RC
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair
select.l<-sum(bNTI.v<(-2))/num.pair
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res,"gen_Processes.csv")

##plot
#Data Loading
bnti<-read.csv("spe.weighted_bNTI.csv",row.names = 1)
process<-read.csv("gen_Processes.csv",row.names = 1)
p1=ggplot(df,aes(value,fill=mp, color=mp)) +
  xlab("bNTI") + ylab("Density")+
  geom_density(alpha = 0.2,size=1.2)+geom_rug()+
  scale_fill_manual(values=c("#B3808B"))+ 
  scale_color_manual(values=c("#B3808B"))+
  theme_bw()
p2=ggplot(process,aes(Process, value),position="stack") +
  geom_bar(aes(fill = Process), stat = "identity",color="black",size=0.4,
           position = "fill", width = 0.6,data=process)+
  scale_fill_manual(values=c("#636270", "#CB4D23","#E9244A","#3C6CAF", "#E33F37"))+
  theme_bw()
