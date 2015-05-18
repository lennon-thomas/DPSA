SummaryPanel<- function(AssessData,LengthDat,Species,Site,YearsToSmooth,Theme)
{
    

  
#     LengthDat<- LengthData
#     AssessData<- AssessmentResults
  LengthDat$MPA<- as.factor(LengthData$MPA)
  
  MaxYear<- max(LengthDat$Year,na.rm=T)
  
  LengthDat$Year<- as.factor(LengthData$Year)
  
  levels(LengthDat$MPA)<- c('Fished','MPA')
  
  MeanLength<- ddply(LengthDat,c('Year','MPA'),summarize,MeanLength=mean(Length,na.rm=T))
  
  #   xyplot(MeanLength~ Year,group= MPA,data=MeanLength,type='b',lwd=3,auto.key=T)
    
#   AssessData$LowerCI<- AssessData$LowerCI - AssessData$Value
# 
#   AssessData$UpperCI<- AssessData$UpperCI - AssessData$Value
#   
#   
  pdf(file=paste(FigureFolder,'Assessment Summary.pdf',sep=''),width=16,height=14)
  
#   theme(legend.position='top')
  LengthPlot<- (ggplot(data=subset(LengthDat,Year==MaxYear),aes(Length,fill=MPA))+
                 geom_density(alpha=0.7,aes(y=..count..))+xlab('Length (cm)')
                +geom_vline(xintercept=Fish$Mat50,linetype='longdash',size=2)
                +scale_fill_manual(name='',values=c(FishedColor,MPAColor))+Theme+ylab('Count'))
  
  limits <- aes(ymax = UpperCI, ymin=LowerCI)
  
  AssessData$LowerCI[is.na(AssessData$LowerCI)]<- AssessData$Value[is.na(AssessData$LowerCI)]

  AssessData$UpperCI[is.na(AssessData$UpperCI)]<- AssessData$Value[is.na(AssessData$UpperCI)]
#   FPlot<- (ggplot(data=subset(AssessData,Method=='CatchCurve'),aes(x=Year,y=Value))+
#     geom_smooth(se=F,size=2)+geom_errorbar(limits)+ylab('F/Fmsy')+geom_hline(yintercept=1))
  FPlot<- (ggplot(data=subset(AssessData,(Method=='CatchCurve' | Method=='LBSPR') & Metric=='FvM'),aes(x=Year,y=Value,color=Method))+
             geom_smooth(se=F,aes(ymin=LowerCI,ymax=UpperCI,fill=Method),stat='identity',size=2,alpha=0.4,size=2)+ylab('F/M')+
             geom_hline(yintercept=1,size=2,linetype='longdash')
           +scale_color_manual(name='',values=c("#1f78b4","#33a02c"))+scale_fill_manual(name='',values=c("#1f78b4","#33a02c"))+Theme)
  
  
  DensityPlot<- (ggplot(data=subset(AssessData,Method=='DensityRatio'),aes(x=Year,y=Value))+
             geom_smooth(se=F,stat='identity',aes(ymin=LowerCI,ymax=UpperCI),fill='#253494',size=2,color='#253494')+ylab('Density Ratio')+geom_hline(yintercept=1,size=2,linetype='longdash')
             +Theme)
  
  SPRPlot<- (ggplot(data=subset(AssessData,Method=='LBSPR' & Metric=='SPR'),aes(x=Year,y=Value))+
                   geom_smooth(se=F,stat='identity',aes(ymin=LowerCI,ymax=UpperCI),size=2,fill='#b30000',color='#b30000')+ylab('SPR')+geom_hline(yintercept=0.4,size=2,linetype='longdash')
             +Theme)
  
grid.arrange(LengthPlot,FPlot,DensityPlot,SPRPlot
    ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
  dev.off()

pdf(file=paste(FigureFolder,'Assessment Summary 2.pdf',sep=''),width=16,height=14)

#   theme(legend.position='top')
LengthPlot<- (ggplot(data=subset(LengthDat,Year==MaxYear),aes(Length,fill=MPA))+
                geom_density(alpha=0.7,aes(y=..count..))+xlab('Length (cm)')
              +geom_vline(xintercept=Fish$Mat50,linetype='longdash',size=2)
              +scale_fill_manual(name='',values=c(FishedColor,MPAColor))+Theme+ylab('Count'))

limits <- aes(ymax = UpperCI, ymin=LowerCI)

AssessData$LowerCI[is.na(AssessData$LowerCI)]<- AssessData$Value[is.na(AssessData$LowerCI)]
AssessData$UpperCI[is.na(AssessData$UpperCI)]<- AssessData$Value[is.na(AssessData$UpperCI)]
#   FPlot<- (ggplot(data=subset(AssessData,Method=='CatchCurve'),aes(x=Year,y=Value))+
#     geom_smooth(se=F,size=2)+geom_errorbar(limits)+ylab('F/Fmsy')+geom_hline(yintercept=1))
FPlot<- (ggplot(data=subset(AssessData,(Method=='CatchCurve') & Metric=='FvM'),aes(x=Year,y=Value))+
           geom_smooth(se=F,aes(ymin=LowerCI,ymax=UpperCI),stat='identity',size=2,alpha=0.4,size=2,
                       color="#336600",fill="#336600")+ylab('F/M')+
           geom_hline(yintercept=1,size=2,linetype='longdash')+Theme)

DensityPlot<- (ggplot(data=subset(AssessData,Method=='DensityRatio'),aes(x=Year,y=Value))+
                 geom_smooth(se=F,stat='identity',aes(ymin=LowerCI,ymax=UpperCI),fill='#253494',size=2,color='#253494')+ylab('Density Ratio')+geom_hline(yintercept=1,size=2,linetype='longdash')
               +Theme)

SPRPlot<- (ggplot(data=subset(AssessData,Method=='LBSPR' & Metric=='SPR'),aes(x=Year,y=Value))+
             geom_smooth(se=F,stat='identity',aes(ymin=LowerCI,ymax=UpperCI),size=2,fill='#b30000',color='#b30000')+ylab('SPR')+geom_hline(yintercept=0.4,size=2,linetype='longdash')
           +Theme)

grid.arrange(LengthPlot,FPlot,DensityPlot,SPRPlot
             ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
dev.off()
  
}