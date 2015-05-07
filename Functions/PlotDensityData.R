PlotDensityData<- function(DensityDat,FigureFolder,Fish,Species,Site)
{

#   DensityDat<- DensityData
    
  DensitySummary<- ddply(DensityDat,c('Year','MPA'),summarize,NumberDensity=mean(Count/SampleArea,na.rm=T),BiomassDensity=mean(Biomass/SampleArea,na.rm=T))
  
  DensitySummary$SiteType[DensitySummary$MPA==0]<- 'Fished'
  
  DensitySummary$SiteType[DensitySummary$MPA==1]<- 'MPA'
  
  pdf(file=paste(FigureFolder,Species,'-',Site,' Density Data.pdf',sep=''))
  print(ggplot(data=DensitySummary,aes(Year,NumberDensity,color=SiteType))+geom_line(size=3)+
          scale_color_manual(name='',values=c(FishedColor,MPAColor))+
    ylab(expression(paste('Number per ',km^{2},sep='')))+ggtitle(paste(Species,Site,sep='-')))
  
  print(ggplot(data=DensitySummary,aes(Year,BiomassDensity,color=SiteType))+geom_line(size=3)+
          scale_color_manual(name='',values=c(FishedColor,MPAColor))+
     ylab(expression(paste('Biomass per ',km^{2},sep='')))+ggtitle(paste(Species,Site,sep='-')))
  
  
  print((ggplot(data=DensityDat,aes(Year,Biomass/SampleArea,color=MPA)))+geom_point(size=3,alpha=0.6)
  +geom_smooth(size=2,se=F)+ylab(expression(paste('Biomass/',km^{2}, ' per trip',sep='')))+ggtitle(paste(Species,Site,sep='-'))
  +scale_color_manual(name='',values=c(FishedColor,MPAColor)))
  
  dev.off()
  write.csv(file=paste(ResultFolder,AssessmentName,' Density Data Summary.csv',sep=''),DensitySummary)
  
}