###### CCFRP Data Poor Stock Assessments ######
#This wrapper calls the AssessmentModules and SubFunctions of the DPSA modules 
#created by Ovando et al. to run DPSAs on the CCFRP groundfish data

#Created by Dan Ovando


### Setup Working Environment ###
rm(list=ls())
set.seed(446)

# library(lattice)
library(plyr)
# library(grid)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(ggmap)
library(animation)
library(R2admb)
library(tidyr)


sapply(list.files(pattern="[.]R$", path="Functions", full.names=TRUE), source)

source("AssessmentModules.R") #Pull in assessment modules
source("SubFunctions.R") #Pull in helper functions for assessment modules

# High Level Assessment Options -------------------------------------------

Assessment<- 'CCFRP 8.2'

NumIterations<- 10

Font<- 'Helvetica'

FontColor<- 'Black'

RunAssessments<- TRUE

NumberOfSpecies<- 5

ReserveYear<- 2007

Assessments<- c('CatchCurve','DensityRatio','LBSPR')

DefaultSD<- 0.001

MinSampleSize<- 150

MPAColor<- "#1b9e77"

FishedColor<- "#d95f02"

### Pull in Assessment Data ###

GFD<- read.csv('CCFRP Data/For_Jono_CCFRPdata_April2014.csv',stringsAsFactors=F) #Read GroundFishData (GFD)

GFD<- subset(GFD,is.na(Year)==F)

Locations<- read.csv('CCFRP Data/CCFRP_2014Data_withGPS.csv')

Locations$SiteId<- paste(Locations$Area,Locations$Site..MPA..REF,sep='-')

SiteGPS<- ddply(Locations,c('SiteId'),summarize,MeanLon=mean(Lon.Center.Point,na.rm=T),MeanLat=mean(Lat.Center.Point,na.rm=T))

GFD$SiteId<- paste(GFD$Site,GFD$MPA_or_REF,sep='-')

GFD<- plyr::join(GFD,SiteGPS,by='SiteId')

SpeciesNames<- read.csv('CCFRP Data/Fish Species.csv',stringsAsFactors=F) #Read species names

LifeHistory<- read.csv('CCFRP Data/CCFRP Life History 2.csv',stringsAsFactors=F) #Read life history

LHI<- read.csv('CCFRP Data/Prince et al 2014 LHI table.csv',stringsAsFactors=F) %>% subset(Taxa=='Teleost ') #Read Prince et al. 2014 LHI

LifeHistory<- LifeHistory[,2:dim(LifeHistory)[2]]

LifeData<- colnames(LifeHistory)[3:dim(LifeHistory)[2]]

source('CCFRP Data/Default_Controlfile.R')

GFD<- join(GFD,SpeciesNames,by='Species.Code')

GFD<- join(GFD,LifeHistory,by='CommName')

GFD<- FindFishbase(GFD)

GFD$AgeMat50<- NA

GFD$AgeMatSource<- NA

GFD$t0[is.na(GFD$t0)]<- 0

SpeciesCatches<- ddply(GFD,c('CommName'),summarize,NumberSampled=length(length_cm),HasLifeHistory=mean(vbk))

SpeciesCatches$NumberSampled<- SpeciesCatches$NumberSampled*as.numeric(is.na(SpeciesCatches$HasLifeHistory)==F)*as.numeric((SpeciesCatches$NumberSampled)>MinSampleSize)

SpeciesCatches<- SpeciesCatches[order(SpeciesCatches$NumberSampled,decreasing=T),]

TopSpecies<- SpeciesCatches$CommName[SpeciesCatches$NumberSampled>0 & is.na(SpeciesCatches$HasLifeHistory)==F]

GFD<- GFD[GFD$CommName %in% TopSpecies,]

GFD<- AddMissingFish(GFD)

# Sites<- c('All',unique(GFD$Site))

Sites<- c('All')

AssessmentResults<- list()

Counter<- 0

if (RunAssessments==T)
{
  
  for (s in 1:length(Sites))    
  {
    
    show(Sites[s])
    WhereSite<- GFD$Site==Sites[s]
    
    if (Sites[s]=='All'){WhereSite<- rep(1,dim(GFD)[1])==1}
    
    
    TopSpecies<- ddply(GFD[WhereSite,],c('CommName'),summarize,NumSamples=length(length_cm[length_cm>0]))
    
    Fishes<- subset(TopSpecies,NumSamples>MinSampleSize)$CommName
    
      for (f in 1:length(Fishes))
#     for (f in 2)
    {
      
      show(Fishes[f])
      
      Species<- Fishes[f]
      
      iGFD<- GFD[WhereSite & GFD$CommName==Fishes[f],]
      
      AssessmentName <- paste(Assessment,Sites[s],Fishes[f],sep='_')
      
      Directory<- paste(Assessment, "/", Sites[s],'/',Fishes[f],'/',sep='')
      
      if (file.exists(Directory)==F)
      {
        dir.create(Assessment)
        dir.create(paste(Assessment, "/", Sites[s],sep=''))
        dir.create( paste(Assessment, "/", Sites[s],'/',Fishes[f],'/',sep=''))
      }
      
      # Format Data -------------------------------------------------------------
      
      SpeciesLifeHistory<- iGFD[1,colnames(iGFD) %in% LifeData]
      
      HasLifeHistory<- SpeciesLifeHistory[which(is.na(SpeciesLifeHistory)==F)]
      
      HasLife<- colnames(HasLifeHistory)
      
      for (l in 1:length(HasLife))
      {
        WhereLife<- which(names(Fish)==HasLife[l])
        
        Fish[[WhereLife]]<- as.numeric(HasLifeHistory[l])  
      }
      
      Fish$M<- Fish$vbk*Fish$MvK
      
      #     if (is.na(Fish$MaxAge))
      #     {
      Fish$MaxAge<- ceiling(-log(0.005)/Fish$M)
      #     }
      #     
      ClosestMvK<- which((Fish$MvK-LHI$MeanMvK)^2==min((Fish$MvK-LHI$MeanMvK)^2))[1]
      
      Fish$LengthMatRatio<- LHI$MeanLMATvMLinf[ClosestMvK]
      
      Fish$MinMvK<- LHI$MinMvK[ClosestMvK]
      
      Fish$MaxMvK<- LHI$MaxMvK[ClosestMvK]
      
      Fish$MinLengthMatRatio<- LHI$MinLMATvMLinf[ClosestMvK]
      
      Fish$MaxLengthMatRatio<- LHI$MaxLMATvMLinf[ClosestMvK]
      
      if (is.na(SpeciesLifeHistory$Mat50) | SpeciesLifeHistory$Mat50>(SpeciesLifeHistory$Linf* SpeciesLifeHistory$LengthMatRatio) ) #Use Prince et al. 2014 LHI
      {
        Fish$Mat50<- Fish$Linf*Fish$LengthMatRatio
        
        Fish$Mat95<- as.numeric(1.1*Fish$Mat50)
      }
      
      if (is.na(SpeciesLifeHistory$AgeMat50)==F)
      {
        Fish$Mat50<- LengthAtAge(SpeciesLifeHistory$AgeMat50,Fish,0)
        
        Fish$Mat95<- as.numeric(1.1*Fish$Mat50)
      }
      
      Fish$Mat95<- as.numeric(1.1*Fish$Mat50)
      
      ReformData<- FormatCCFRPData(iGFD)
      
      LengthData<- ReformData$LengthData
      
      DensityData<- ReformData$DensityData
      
      write.csv(file=paste(Directory,AssessmentName,'_LengthData.csv',sep=''),LengthData)
      
      write.csv(file=paste(Directory,AssessmentName,'_DensityData.csv',sep=''),DensityData)
      
      FigureFolder<- paste(Directory,'Figures/',sep='')
      
      ResultFolder<- paste(Directory,'Results/',sep='')
      
      if (file.exists(FigureFolder)==F)
      {
        dir.create(FigureFolder,recursive=T)
        
        dir.create(ResultFolder,recursive=T)
      }
      
      write.csv(file=paste(ResultFolder,Species,' Life History.csv',sep=''),as.data.frame(Fish))
      
      PlotLifeHistory()
      
      Theme<- theme(legend.position='top',plot.background=element_rect(color=NA),
                    rect=element_rect(fill='transparent',color=NA)
                    ,text=element_text(size=40,family=Font,color=FontColor),
                    axis.text.x=element_text(color=FontColor),
                    axis.text.y=element_text(color=FontColor),legend.key.size=unit(2.5,'cm'))
      
      PlotLengthData(LengthData,FigureFolder,Fish,Species,Sites[s],Theme)
      
      PlotDensityData(DensityData,FigureFolder,Fish,Species,Sites[s],Theme)
      
      #     MapCCFRP(ReformData)
      
      ### Run Assessments ###
      
      #     AssessmentResults<- as.data.frame(matrix(NA,nrow=length(Assessments)*10,ncol=9))
      #     
      #     colnames(AssessmentResults)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
      #     
      #     AssessmentResults$Year<- as.numeric(AssessmentResults$Year)
      
      Fish$LHITol<- 0.99
      
      for (a in 1:length(Assessments)) #Loop over possible assessments, store in Assessment results. Many assessments have more detailed outputs than can also be accessed 
      {
        
        Counter<- Counter+1
        if (Assessments[a]=='LBAR') #Run LBAR assessment
        {
          
          SampleCheck<- CheckLengthSampleSize(LengthData)        
          
          if (SampleCheck$YearsWithEnoughData>0)
          {
            
            
            Temp<- LBAR(SampleCheck$ParedData,LagLength=1,Weight=0.2,IncludeMPA=0,ReserveYr=ReserveYear,OutsideBoundYr=NA,Iterations=NumIterations,
                        BootStrap=1,LifeError=1,Lc=NA)$Output		
            
            StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
              rename(Site=Sites.s.,Assessment=Assessments.a.)
            
            AssessmentResults[[Counter]]<-StoreAssess
          }
        }
        
        if (Assessments[a]=='CatchCurve') #Run Catch Curve analysis
        {
          
          SampleCheck<- CheckLengthSampleSize(LengthData)        
          
          if (SampleCheck$YearsWithEnoughData>0)
          {
            
            
            Temp<- CatchCurve(SampleCheck$ParedData,CatchCurveWeight='AgeBased',WeightedRegression=1,
                              ReserveYr=ReserveYear,OutsideBoundYr=NA,ManualM=0,Iterations=NumIterations,BootStrap=1,LifeError=1,HistInterval=1)$Output
            
            StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
              rename(Site=Sites.s.,Assessment=Assessments.a.)
            
            AssessmentResults[[Counter]]<-StoreAssess
            
          }
        }
        
        
        if (Assessments[a]=='DensityRatio') #Run density ratio analysis 
        {
          
          Temp<- DensityRatio(DensityData,LagLength=1,Weight=1,Form='Biomass',Iterations=NumIterations,BootStrap=1)$Output
          
          StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
            rename(Site=Sites.s.,Assessment=Assessments.a.)
          
          AssessmentResults[[Counter]]<-StoreAssess
        }
        
        if (Assessments[a]=='CatchMSY')
        {
          Temp2<- CatchMSY(CatchData,1000,0.05,0,0,1,0,0,1,NA,c(0.75,0.99),NA,NA,c(0.25,0.65))
          
          Temp<- Temp2$Output
          
          StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
            rename(Site=Sites.s.,Assessment=Assessments.a.)
          
          AssessmentResults[[Counter]]<-StoreAssess
          
        }
        
        if (Assessments[a]=='LBSPR') #Run LBSPR Assessment
        {
          
          SampleCheck<- CheckLengthSampleSize(LengthData)        
          
          if (SampleCheck$YearsWithEnoughData>0)
          {
            
            LengthQuantile<- quantile(SampleCheck$ParedData$Length,na.rm=T)
            
            
            #           Temp2<- LBSPR(SampleCheck$ParedData,EstimateM=0,Iterations=1,BootStrap=1,
            #                         LifeError=1,LengthBins=1,ReserveYear=ReserveYear,SL50Min=LengthQuantile[1],
            #                         SL50Max=LengthQuantile[2],DeltaMin=NA,DeltaMax=NA,IncludeReserve=TRUE)
            
            Temp2<- LBSPR(SampleCheck$ParedData,EstimateM=0,Iterations=NumIterations,BootStrap=1,
                          LifeError=1,LengthBins=1,ReserveYear=ReserveYear,SL50Min=LengthQuantile[1],
                          SL50Max=LengthQuantile[2],DeltaMin=0.01,DeltaMax=.2*Fish$Linf,IncludeReserve=FALSE)
            
            
            Temp<- Temp2$Output
            
            StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
              rename(Site=Sites.s.,Assessment=Assessments.a.)
            
            AssessmentResults[[Counter]]<-StoreAssess
          }
        }
        
        
        if (Assessments[a]=='DBSRA') #Run DBSRA Assessment
        {
          
          DCAC.start.yr <- CatchData$Year[1] #start of the catch period
          DCAC.end.yr<- CatchData$Year[length(CatchData$Year)] #end of the catch period
          delta.yr<- CatchData$Year[length(CatchData$Year)] #Year that current depletion is fit to
          DBSRA.OFL.yr<- CatchData$Year[length(CatchData$Year)] #Year to calculate DBSRA OFL outputs
          FMSYtoMratio <- 0.8 #ratio of Fmsy to M
          SD.FMSYtoMratio<- 0.05
          Delta<- 0.7
          SD.Delta<- 0.1
          DeltaLowerBound<- 0.5
          DeltaUpperBound<- 0.9
          BMSYtoB0ratio <- 0.3
          SD.BMSYtoB0ratio<- 0.1
          BMSYtoB0LowerBound<- 0.2
          BMSYtoB0UpperBound<- 0.5
          CatchInterp<-1
          NIter<- 500	
          
          
          Temp2<- DBSRA(CatchData, DCAC.start.yr, DCAC.end.yr, delta.yr, DBSRA.OFL.yr, FMSYtoMratio, SD.FMSYtoMratio, Delta, SD.Delta, DeltaLowerBound, DeltaUpperBound, BMSYtoB0ratio, SD.BMSYtoB0ratio, BMSYtoB0LowerBound, BMSYtoB0UpperBound, CatchInterp, NIter)
          
          Temp<- Temp2$Output
          
          StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
            rename(Site=Sites.s.,Assessment=Assessments.a.)
          
          AssessmentResults[[Counter]]<-StoreAssess
        }
        
        
        show(paste('Finished ',Assessments[a],'-',round(100*a/length(Assessments),2),'% Done',sep=''))
        
      }
      
      CurrentResults<- ldply(AssessmentResults) %>% subset(Species==Fishes[f] & Site==Sites[s])
      
      if (any(CurrentResults$Assessment=='CatchCurve') & any(CurrentResults$Assessment=='DensityRatio') 
          & any(CurrentResults$Assessment=='LBSPR') )
      {
        
        Theme<- theme(legend.position='top',plot.background=element_rect(color=NA),
                      rect=element_rect(fill='transparent',color=NA)
                      ,text=element_text(size=30,family=Font,color=FontColor),
                      axis.text.x=element_text(color=FontColor),
                      axis.text.y=element_text(color=FontColor))
        
        SummaryPanel(CurrentResults,LengthData,Species,Sites[s],YearsToSmooth=3,Theme)
        
      }
      save.image(file=paste(ResultFolder,AssessmentName,'_Settings.RData',sep='')) #Save settings used to produce current results
      write.csv(file=paste(ResultFolder,AssessmentName,'_Results.csv',sep=''),CurrentResults) #Save current results
      
    } #Close species  (f) loop
  } #Close sites (s) loop
  
  
  save(file=paste(Assessment,'/Assessment Results.Rdata',sep=''),AssessmentResults)
  
  save(file=paste(Assessment,'/Raw Data.Rdata',sep=''),GFD)
}
if (RunAssessments==F)
{
  try(load(file=paste(Assessment,'/Assessment Results.Rdata',sep='')))
  try(load(file=paste(Assessment,'/Raw Data.Rdata',sep='')))
}

FlatAssessments<- ldply(AssessmentResults)


SAData<- read.csv('CCFRP Data/CCFRP Stock Assessment Values.csv')

SAData$Assessment<- 'Stock Assessment'

ddply(FlatAssessments,c('Species','Year','Metric'),summarize,Number=mean(Value))



Ramp<- 'RdYlGn'

LineKeynoteTheme<- theme(plot.background=element_rect(color=NA),rect=element_rect(fill='transparent',color=NA)
                         ,text=element_text(size=30,family=Font,color=FontColor))


KeynoteTheme<- theme(plot.background=element_rect(color=NA),rect=element_rect(fill='transparent',color=NA)
                     ,text=element_text(size=30,family=Font,color=FontColor),
                     axis.text.x=element_text(angle=45,hjust=0.9,vjust=0.9,color=FontColor),
                     axis.text.y=element_text(color=FontColor),legend.position='none')

# TimeToNPB$Species <- reorder(TimeToNPB$Species, 1/TimeToNPB$BenefitYear)

AggPlot<- function(PlotDat,PlotName,Theme,Order)
{
  
  PlotDat$Species<- as.factor(PlotDat$Species)
  
  PlotDat$Species<- reorder((PlotDat$Species),PlotDat$Value)
  
  if (Order=='Descending')
  {
    PlotDat$Species<- reorder((PlotDat$Species),1/(PlotDat$Value+abs(min(PlotDat$Value,na.rm=T))))
  }
    
    

    print(ggplot(data=PlotDat,
                 aes(x=(Species),y=Value,fill=Species))+scale_fill_brewer(palette=Ramp)
          +ylab(PlotName)+xlab(NULL)
          +geom_bar(stat='identity',color='black')+
            Theme+geom_hline(yintercept=0,size=2))

}

OrderDat<- function(PlotDat,OrdVar,OrdBy,Order)
{
  PlotDat[,OrdVar]<- as.factor(PlotDat[,OrdVar])
  
  PlotDat[,OrdVar]<- reorder(PlotDat[,OrdVar],PlotDat[,OrdBy])
  
  if (Order=='Descending')
  {
    PlotDat$Species<- reorder(PlotDat[,OrdVar],1/(PlotDat[,OrdBy]+abs(min(PlotDat[,OrdBy],na.rm=T))))
  }
  
  return(PlotDat)
}

pdf(paste(Assessment,'/Across Species Summary.pdf',sep=''),width=18,height=12)

PlotDat<- subset(FlatAssessments,Year==2014 & Metric=='FvM' & Method=='CatchCurve' & Site=='All')

AggPlot(PlotDat,'FvM',KeynoteTheme,'Descending')

PlotDat<- subset(FlatAssessments,Year==2014 & Metric=='FvM' & Method=='LBSPR' & Site=='All')

AggPlot(PlotDat,'FvM',KeynoteTheme,'Descending')

PlotDat<- subset(FlatAssessments,Year==2014 & Metric=='SPR' & Site=='All')

AggPlot(PlotDat,'SPR',KeynoteTheme,'Ascending')

PlotDat<- subset(FlatAssessments,Year==2014 & Metric=='Biomass Density Ratio' & Site=='All')

AggPlot(PlotDat,'Density Ratio',KeynoteTheme,'Ascending')

AllAssess<- subset(FlatAssessments,Site=='All')  %>% mutate(SA=paste(Species,Method,Metric,sep='-'))

SAS<- unique(AllAssess$SA)

AllAssess$FSlope<- NA

for (i in 1:length(SAS))
{
  
  FSlope<- lm(Value~Year,data=subset(AllAssess,SA==SAS[i]))$coefficient[2]
  
  AllAssess$Slope[AllAssess$SA==SAS[i]]<- FSlope
}

Ordered<- OrderDat(subset(AllAssess,Year==2014 & Metric=='FvM' & Method=='CatchCurve'),'Species','Slope','Descending')

print(ggplot(data=,Ordered,aes(x=factor(Species),y=Slope,fill=Species))+scale_fill_brewer(palette=Ramp)+ylab('Trend in F/M')
      +geom_bar(stat='identity',color='black')+KeynoteTheme+xlab(NULL)+geom_hline(yintercept=0,size=2))


Ordered<- OrderDat(subset(AllAssess,Year==2014 & Metric=='SPR'),'Species','Slope','Increasing')


print(ggplot(data=,Ordered,aes(x=factor(Species),y=Slope,fill=Species))+scale_fill_brewer(palette=Ramp)+ylab('Trend in SPR')
      +geom_bar(stat='identity',color='black')+KeynoteTheme+xlab(NULL)+geom_hline(yintercept=0,size=2))


Ordered<- OrderDat(subset(AllAssess,Year==2014 & Metric=='Biomass Density Ratio'),'Species','Slope','Increasing')

print(ggplot(data=,Ordered,aes(x=factor(Species),y=Slope,fill=Species))+scale_fill_brewer(palette=Ramp)+ylab('Trend in DR')
      +geom_bar(stat='identity',color='black')+KeynoteTheme+xlab(NULL)+geom_hline(yintercept=0,size=2))

dev.off()

