PlotLengthData<- function(LengthDat,FigureFolder,Fish,Species,Site)
{
#   LengthDat<- LengthData
  
  LengthDat<- subset(LengthDat,is.na(Length)==F & Length>0)
  
  LengthDat$SiteType[LengthDat$MPA==0]<- 'Fished'

  LengthDat$SiteType[LengthDat$MPA==1]<- 'MPA'
  
  pdf(file=paste(FigureFolder,Species,'-',Site,'Length Data.pdf',sep=''))
  print(ggplot(data=LengthDat,aes(Length,fill=(SiteType)))+geom_density(alpha=0.6,aes(y=..count..))+scale_fill_discrete(name='')
  +geom_vline(xintercept=Fish$Mat50,linetype='longdash')+
    geom_vline(xintercept=Fish$Linf,linetype='longdash',color='red2')+
facet_wrap(~Year,as.table=F)+ggtitle(paste(Species,Site,sep='-')))
  
print(ggplot(data=LengthDat,aes(factor(Year),Length,fill=SiteType))+geom_boxplot(varwidth=F,alpha=1))
  dev.off()
  
}