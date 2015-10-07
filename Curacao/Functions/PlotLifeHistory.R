PlotLifeHistory<- function()
{
  
  Lengths<- LengthAtAge(1:Fish$MaxAge,Fish,0)
  
  Maturity<- MaturityAtAge(Lengths,Fish)
  
  LifeMat<- data.frame(1:Fish$MaxAge,Lengths,Maturity)
  
  colnames(LifeMat)<- c('Age','Length','PercentMature')
  
  pdf(file=paste(FigureFolder,AssessmentName,' Life History.pdf',sep=''))
  
  LifePlot<- ggplot(LifeMat,aes(x=Age,y=Length,color=PercentMature))+geom_line(size=2)+ylab('Length (cm)')
  LifePlot<-LifePlot+theme(axis.line=element_line(size=0.6,colour="black"),panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank())
  
  print( LifePlot + scale_colour_gradient(limits=c(0, 1), low="steelblue2",high='red'))  
  
  dev.off()
}
