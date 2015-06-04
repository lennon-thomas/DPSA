PlotCPUEData<- function(CPUEDat,FigureFolder,Fish,Species,Site,Theme)
{
  
  #   DensityDat<- DensityData
  CPUESummary<- ddply(CPUEDat,c('Year','MPA'),summarize,NumberCPUE=mean(Count/AnglerHours,na.rm=T),BiomassCPUE=mean(Biomass/AnglerHours,na.rm=T))
  
  CPUESummary$SiteType[CPUESummary$MPA==0]<- 'Fished'
  
  CPUESummary$SiteType[CPUESummary$MPA==1]<- 'MPA'
  
  CPUEDat$SiteType[CPUESummary$MPA==0]<- 'Fished'
  
  CPUEDat$SiteType[CPUESummary$MPA==1]<- 'MPA'
  
  pdf(file=paste(FigureFolder,Species,'-',Site,' CPUE Data.pdf',sep=''),height=5,width=5)
  print(ggplot(data=CPUESummary,aes(Year,NumberCPUE,color=SiteType))+geom_line(size=3)+
          scale_color_manual(name='',values=c(FishedColor,MPAColor))+
          ylab(paste('Number/Angler Hour',sep=''))+ggtitle(paste(Species,Site,sep='-'))+Theme)
  
  print(ggplot(data=CPUESummary,aes(Year,BiomassCPUE,color=SiteType))+geom_line(size=3)+
          scale_color_manual(name='',values=c(FishedColor,MPAColor))+
          ylab(paste('kg/Angler Hour',sep=''))+ggtitle(paste(Species,Site,sep='-'))+Theme)
  
  print((ggplot(data=CPUEDat,aes(Year,Biomass/AnglerHours,color=SiteType)))+geom_point(size=6,alpha=0.6)
        +ylab(paste('Number/Angler Hour per trip',sep=''))+ggtitle(paste(Species,Site,sep='-'))
        +scale_color_manual(name='',values=c(FishedColor,MPAColor))+Theme)
  
  dev.off()
  write.csv(file=paste(ResultFolder,AssessmentName,' CPUE Data Summary.csv',sep=''),CPUESummary)
  
}