library(directlabels)
data<-read.csv("results3.csv")
names(data)
data$Species[1]

p<-ggplot(data=data,aes(x=bio,y=fish))+
  geom_hline(yintercept=mean(1),linetype="dashed",size=0.6)+
  geom_vline(xintercept=mean(1),linetype="dashed",size=0.6)+
  geom_point(aes(colour=Species,shape=Assessment.Method),size=c(3))+
  labs(x=expression(B/B[MSY]),y=expression(F/F[MSY]))+
  scale_shape_manual(values=c(16,10))+
  scale_y_continuous(limits=c(0,7), expand = c(0, 0))+scale_x_continuous(limits=c(0,2.5))
p<-p+theme(
  axis.text = element_text(size = 14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white"),
  axis.line=element_line(colour="Black",size=0.85),axis.text.x = element_text(size=14),axis.title.y=element_text(size=16,vjust=1.5,face="italic"),
  axis.title.x=element_text(size=16,vjust=0.1,face="italic"),legend.position="right",legend.title=element_text(size=14),legend.text=element_text(size=12))
p<-p+coord_fixed(.4)

plot(p)

pdf("kobeplot.pdf")
p
dev.off()
