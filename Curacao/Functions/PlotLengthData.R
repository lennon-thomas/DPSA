PlotLengthData<- function(LengthDat,FigureFolder,Fish,Species,Site,Theme)
{
#    LengthDat<- LengthData
#   Site<-"Curacao"
#  Species<-"Acanthurus chirurgus"
  LengthDat<- subset(LengthDat,is.na(Length)==F & Length>0)
  
  LengthDat$Year[LengthDat$Year==2008]<- '2008'

  LengthDat$Year[LengthDat$Year==2014]<- '2014'
  
  LengthDat$Year[LengthDat$Year==2015]<- '2015'
  
  pdf(file=paste(FigureFolder,Species,'-',Site,'Length Data.pdf',sep=''),width=7,height=5)
  
  LengthDist<- (ggplot(data=LengthDat,aes(Length,fill=(Year)))+geom_density(alpha=0.6) # aes(y=..count..))
  +scale_fill_manual(name='',values=c("lightseagreen","orange","blue"))
  +geom_vline(xintercept=Mat50,linetype='longdash')+
   geom_vline(xintercept=Linf,linetype='longdash',color='red2')
   +ggtitle(paste(Species,sep=""))+Theme+xlab("Length(cm)")+ylab("Frequency"))
  
  LengthDist<-LengthDist+scale_fill_discrete(name="Sample Size\n",breaks=c("2008","2014","2015"),labels=c("\nn=134\n ","\nn=398\n ","\nn=28\n "))
  
  #+xlab("Length (cm))"+ylab("Density")
  #facet_wrap(~Year,as.table=F)
  
  print(LengthDist)
  
# pdf(file=paste(FigureFolder,Species,'-',Site,'Length Data.pdf',sep=''))
# print(ggplot(data=LengthDat,aes(Length,fill=(SiteType)))+geom_density(alpha=0.6,aes(y=..count..))
#       +scale_fill_manual(name='',values=c(FishedColor,MPAColor))
#       +geom_vline(xintercept=Fish$Mat50,linetype='longdash')+
#         geom_vline(xintercept=Fish$Linf,linetype='longdash',color='red2')+
#         facet_wrap(~Year,as.table=F)+ggtitle(paste(Species)))
 dev.off()

LengthBox<- (ggplot(data=LengthDat,aes(factor(Year),Length,fill=Year))
+geom_boxplot(varwidth=F,alpha=1)+scale_fill_manual(name='',values=c("lightseagreen","orange","blue"))+Theme+xlab('Year')+
  ylab('Length (cm)'))
  
print(LengthBox)
  dev.off()

save(LengthDist,LengthBox,file=paste(Species,'Length Data.Rdata',sep=''))
    
}