###figure 1A the relative abundance of main phyla
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package
library(dplyr)
library(readxl)
library(gggenomes)
dflink<-read_excel("phylum.xlsx",sheet = "Sheet2")
dflink
gggenomes(genes = df,links = dflink)+
  geom_link(offset = 0.05,aes(fill=group2))+
  geom_gene(shape=0,aes(fill=group),size=5)+
  theme(legend.position = "none")

df3<-read_excel("gggenomes_examples.xlsx",
                sheet = "Sheet3")
df3
df3 %>% 
  filter(seq_id=="MPs") %>% 
  mutate(length=end-start+1) -> dflink01

df3 %>% 
  filter(seq_id=="Soil") %>% 
  mutate(length=end-start+1) -> dflink02

bind_cols(dflink01,dflink02) -> dflinks
colnames(dflinks)<-c(c("seq_id","start","end","group",'length'),
                     paste0(c("seq_id","start","end","group",'length'),2))
dflinks
gggenomes(genes=df3,links = dflinks)+
  geom_bin_label(size=5)+
  geom_link(aes(fill=group),offset = 0.3)+
  geom_gene(aes(fill=group),shape=0,
            size=10)+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank())+
  xlim(-15,100)+
  scale_fill_manual(values =c("#E5A429", "#159ACC", "#0A9C78","#EFE454","#159ACC", "#D55E1D","#CF88AA" ,"#E33F37",
                              "#5FB5E6", "#FA8181", "#5A845B", "#CFCFAA","#0A51A3","#4F4F4F","#663865","#ADC9CA"))+
  annotate(geom = "segment",
           x=0,xend = 100,y = 0.5,yend=0.5)+
  geom_segment(data=data.frame(x=c(0,50,100),
                               xend=c(0,50,100),
                               y=c(0.5,0.5,0.5),
                               yend=c(0.4,0.4,0.4)),
               aes(x=x,xend=xend,y=y,yend=yend),
               inherit.aes = FALSE)
