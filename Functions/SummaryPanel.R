SummaryPanel<- function(AssessData,LengthDat,Species,Site)
{
  
#   LengthDat<- LengthData
#   AssessData<- AssessmentResults
    
  LengthDat$MPA<- as.factor(LengthData$MPA)
  
  MaxYear<- max(LengthDat$Year,na.rm=T)
  
  LengthDat$Year<- as.factor(LengthData$Year)
  
  levels(LengthDat$MPA)<- c('Fished','MPA')
  
  MeanLength<- ddply(LengthDat,c('Year','MPA'),summarize,MeanLength=mean(Length,na.rm=T))
  
#   xyplot(MeanLength~ Year,group= MPA,data=MeanLength,type='b',lwd=3,auto.key=T)
  
pdf(file=paste(FigureFolder,'Assessment Summary.pdf',sep=''))
  grid.arrange(
      densityplot(~Length | Year,groups=MPA,data=LengthDat[LengthDat$Year==MaxYear,],auto.key=T,type='count',lwd=3,panel=function(x,...)
      {
        panel.densityplot(x,...)
        panel.abline(v=Fish$Mat50,lty=2)
      } 
      ),
    
    xyplot(Value ~ Year | Method,group=Metric,data=AssessData,subset=Method=='CatchCurve',type='b',lwd=3,ylab='F/M',
           panel=function(x,y,...)
             {
             panel.xyplot(x,y,...)
             panel.abline(h=1,lty=2)
           })
    ,
    xyplot(Value ~ Year | Method,group=Metric,data=AssessData,subset=Method=='DensityRatio',col='tomato3',type='b',lwd=3,ylab='Fished/MPA',
           panel=function(x,y,...)
           {
             panel.xyplot(x,y,...)
             panel.abline(h=1,lty=2)
           })
    ,
    xyplot(Value ~ Year | Method,group=Metric,data=AssessData,subset=Method=='LBSPR',type='b',lwd=3,ylab='SPR', panel=function(x,y,...)
    {
      panel.xyplot(x,y,...)
      panel.abline(h=0.4,lty=2)
    })
    
    ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
    dev.off()

pdf(file=paste(FigureFolder,'Assessment Summary V2.pdf',sep=''))
grid.arrange(
  densityplot(~Length,groups=MPA,data=LengthDat[LengthDat$Year==MaxYear,],auto.key=T,type='count',lwd=3,panel=function(x,...)
  {
    panel.densityplot(x,...)
    panel.abline(v=Fish$Mat50,lty=2)
  } 
  ),
  
  xyplot(Value ~ Year,data=AssessData,subset=Method=='CatchCurve',type='b',lwd=3,ylab='F/M',xlab='',
         panel=function(x,y,...)
         {
           panel.xyplot(x,y,...)
           panel.abline(h=1,lty=2)
         })
  ,
  xyplot(Value ~ Year,data=AssessData,subset=Method=='DensityRatio',col='tomato3',type='b',lwd=3,ylab='Fished/MPA',xlab='',
         panel=function(x,y,...)
         {
           panel.xyplot(x,y,...)
           panel.abline(h=1,lty=2)
         })
  ,
  xyplot(Value ~ Year,data=AssessData,subset=Method=='LBSPR',type='b',col='forestgreen',lwd=3,ylab='SPR', panel=function(x,y,...)
  {
    panel.xyplot(x,y,...)
    panel.abline(h=0.4,lty=2)
  })
  
  ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
dev.off()

}