###figure 1B the stemp analysis of the phyla in soil and plastisphere
data1 <- read.table("phylum.txt",header = 1,check.names = F,sep = "\t")
data2 <- read.table("analysis.txt",header = 1,check.names = F,sep = "\t")

#Data Processing
df1 <- melt(data1,id.vars = c("sample","group"))
df1$sample <- factor(df1$sample,levels = rev(data1$sample))
df1$group <- factor(df1$group,levels = c("MPs","Soil"))
df1$facet <- rep("Phyla",times=99)
df2 <- melt(data2,id.vars = "sample")
df2$sample <- factor(df2$sample,levels = rev(data1$sample))
df2$group <- rep(c("MPs","Soil"),each=99)

###the first column diagram
p1 <- ggplot(df1,aes(sample,value))+
  geom_bar(aes(fill=group,color=group))#绘制箱线图
coord_flip()+
  scale_fill_manual(values = c("#B3808B","#537874"))+
  scale_color_manual(values = c("#B3808B","#537874"))+
  #Y-axis range
  scale_y_continuous(limits = c(0, 40))+
  #theme
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black",size=10),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=10))+
  #caption
  labs(y="Proportions (%)",x=NULL)+
  #Add group rectangle
  annotate("rect", xmin = 0, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#EEEEEE") +
  annotate("rect", xmin = 6.5, xmax = 20.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#FFFFFF") 

###the second diagram
library(Rmisc) # Ryan Miscellaneous
data3 <- read.table("analysis.txt",header = 1,check.names = F,sep = "\t")

#plot
p2 <- ggplot(data3, aes(sample,value, color = group)) + 
  geom_errorbar(aes(ymin = value- se, ymax = value + se), 
                width = 0,position = position_dodge(0.8),linewidth=0.5) + 
  geom_point(position = position_dodge(0.8),shape=18,size=3)+
  #Shift the X-axis and Y-axis positions
  coord_flip()+
  #theme
  theme_bw()+
  theme(legend.position = c(0.8,0.95),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size=10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=10))+
  #caption
  labs(y="Functional dispersion values",x=NULL,color=NULL)+
  scale_color_manual(values = c("#B3808B","#537874"))+
  #Add group rectangle
  annotate("rect", xmin = 0, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#EEEEEE") +
  annotate("rect", xmin = 6.5, xmax = 20.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#FFFFFF")+
  #Add a significance marker
  annotate('text', label = '**', x =11, y =80, angle=-90, size =5,color="black")+
  geom_segment(x = 10.5, xend = 11.5, y = 78, yend = 78,color = "black", size = 0.8)+
  facet_grid(~ facet)
p2
#Combine
library(aplot) # Decorate a 'ggplot' with Associated Information
p1%>%insert_right(p2,width = 1)


