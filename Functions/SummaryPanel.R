SummaryPanel<- function(AssessData,LengthDat,Species,Site,YearsToSmooth)
{
  
    LengthDat<- LengthData
    AssessData<- AssessmentResults
  
  LengthDat$MPA<- as.factor(LengthData$MPA)
  
  MaxYear<- max(LengthDat$Year,na.rm=T)
  
  LengthDat$Year<- as.factor(LengthData$Year)
  
  levels(LengthDat$MPA)<- c('Fished','MPA')
  
  MeanLength<- ddply(LengthDat,c('Year','MPA'),summarize,MeanLength=mean(Length,na.rm=T))
  
  #   xyplot(MeanLength~ Year,group= MPA,data=MeanLength,type='b',lwd=3,auto.key=T)
    
  pdf(file=paste(FigureFolder,'Assessment Summary.pdf',sep=''))
  
  LengthPlot<- (ggplot(data=subset(LengthDat,Year==MaxYear),aes(Length,fill=MPA))+
                 geom_density(alpha=0.5)+xlab('Length (cm)')+geom_vline(xintercept=Fish$Mat50,color='red',
                                                                        linetype='longdash',size=1.5))
  
  limits <- aes(ymax = UpperCI, ymin=LowerCI)
  
  AssessData$LowerCI[is.na(AssessData$LowerCI)]<- AssessData$Value[is.na(AssessData$LowerCI)]

  AssessData$UpperCI[is.na(AssessData$UpperCI)]<- AssessData$Value[is.na(AssessData$UpperCI)]
  
  
  FPlot<- (ggplot(data=subset(AssessData,Method=='CatchCurve'),aes(x=Year,y=Value))+
    geom_smooth(se=F,size=2)+geom_errorbar(limits)+ylab('F/Fmsy')+geom_hline(yintercept=1))
    
  DensityPlot<- (ggplot(data=subset(AssessData,Method=='DensityRatio'),aes(x=Year,y=Value))+
             geom_smooth(se=F,size=2,color='green4')+geom_errorbar(limits)+ylab('Density Ratio')+geom_hline(yintercept=1))
  
  SPRPlot<- (ggplot(data=subset(AssessData,Method=='LBSPR'),aes(x=Year,y=Value))+
                   geom_smooth(se=F,size=2,color='red2')+geom_errorbar(limits)+ylab('SPR')+geom_hline(yintercept=0.4))
  
grid.arrange(LengthPlot,FPlot,DensityPlot,SPRPlot
    ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
  dev.off()
  
#   pdf(file=paste(FigureFolder,'Assessment Summary V2.pdf',sep=''))
#   grid.arrange(
#     densityplot(~Length,groups=MPA,data=LengthDat[LengthDat$Year==MaxYear,],auto.key=T,type='count',lwd=3,panel=function(x,...)
#     {
#       panel.densityplot(x,...)
#       panel.abline(v=Fish$Mat50,lty=2)
#     } 
#     ),
#     
#     xyplot(movingAverage(Value,n=YearsToSmooth,centered=TRUE) ~ Year,data=AssessData,subset=Method=='CatchCurve',type='b',lwd=3,ylab='F/M',xlab='',
#            panel=function(x,y,...)
#            {
#              panel.xyplot(x,y,...)
#              panel.abline(h=1,lty=2)
#            })
#     ,
#     xyplot(movingAverage(Value,n=YearsToSmooth,centered=TRUE) ~ Year,data=AssessData,subset=Method=='DensityRatio',col='tomato3',type='b',lwd=3,ylab='Fished/MPA',xlab='',
#            panel=function(x,y,...)
#            {
#              panel.xyplot(x,y,...)
#              panel.abline(h=1,lty=2)
#            })
#     ,
#     xyplot(movingAverage(Value,n=YearsToSmooth,centered=TRUE) ~ Year,data=AssessData,subset=Method=='LBSPR',type='b',col='forestgreen',lwd=3,ylab='SPR', panel=function(x,y,...)
#     {
#       panel.xyplot(x,y,...)
#       panel.abline(h=0.4,lty=2)
#     })
#     
#     ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
#   dev.off()
  
}