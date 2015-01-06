###### CCFRP Data Poor Stock Assessments ######
#This wrapper calls the AssessmentModules and SubFunctions of the DPSA modules 
#created by Ovando et al. to run DPSAs on the CCFRP groundfish data

#Created by Dan Ovando


### Setup Working Environment ###
rm(list=ls())
library(lattice)
library(plyr)
library(grid)
library(gridExtra)
library(ggplot2)

sapply(list.files(pattern="[.]R$", path="Functions", full.names=TRUE), source)

source("AssessmentModules.R") #Pull in assessment modules
source("SubFunctions.R") #Pull in helper functions for assessment modules

# High Level Assessment Options -------------------------------------------

Assessment<- 'CCFRP 2014 Scratch 2'

NumberOfSpecies<- 6

ReserveYear<- 2007

AvailableData<- c('LengthData','DensityData')

Assessments<- c('CatchCurve','DensityRatio','LBSPR')

DefaultSD<- 0.001

MinSampleSize<- 200

### Pull in Assessment Data ###


GFD<- read.csv('/Users/danovando/Desktop/Bren/SFG Work/Consulting/TNC/CCFRP/For_Jono_CCFRPdata_April2014.csv',stringsAsFactors=F) #Read GroundFishData (GFD)

SpeciesNames<- read.csv('/Users/danovando/Desktop/Bren/SFG Work/Consulting/TNC/CCFRP/Fish Species.csv',stringsAsFactors=F) #Read species names

LifeHistory<- read.csv('/Users/danovando/Desktop/Bren/SFG Work/Consulting/TNC/CCFRP/CCFRP Life History 2.csv',stringsAsFactors=F) #Read life history

LifeHistory<- LifeHistory[,2:dim(LifeHistory)[2]]

LifeData<- colnames(LifeHistory)[3:dim(LifeHistory)[2]]

source('/Users/danovando/Desktop/Bren/SFG Work/Consulting/TNC/CCFRP/Default_Controlfile.R')

GFD<- join(GFD,SpeciesNames,by='Species.Code')

GFD<- join(GFD,LifeHistory,by='CommName')

GFD<- FindFishbase(GFD)

GFD$AgeMat50<- NA

GFD$AgeMatSource<- NA

GFD$t0[is.na(GFD$t0)]<- 0

SpeciesCatches<- ddply(GFD,c('CommName'),summarize,NumberSampled=length(length_cm),HasLifeHistory=mean(vbk))

SpeciesCatches$NumberSampled<- SpeciesCatches$NumberSampled*as.numeric(is.na(SpeciesCatches$HasLifeHistory)==F)*as.numeric((SpeciesCatches$NumberSampled)>MinSampleSize)

SpeciesCatches<- SpeciesCatches[order(SpeciesCatches$NumberSampled,decreasing=T),]

TopSpecies<- SpeciesCatches$CommName[SpeciesCatches$NumberSampled>0]

# TopSpecies<- SpeciesCatches$CommName[1:NumberOfSpecies]

#  TopSpecies<- LifeHistory$CommName[LifeHistory$HasLifeHistory==1]

GFD<- GFD[GFD$CommName %in% TopSpecies,]

GFD<- AddMissingFish(GFD)

Sites<- c('All',unique(GFD$Site))

for (s in 1:length(Sites))
{
  
  show(Sites[s])
  WhereSite<- GFD$Site==Sites[s]
  
  if (Sites[s]=='All'){WhereSite<- rep(1,dim(GFD)[1])==1}
  
  Fishes<- unique(GFD$CommName[WhereSite])
  
  for (f in 1:length(Fishes))
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
    
    for (l in 1:length(LifeData))
    {
      WhereLife<- which(names(Fish)==LifeData[l])
      
      Fish[[WhereLife]]<- as.numeric(SpeciesLifeHistory[l])  
    }
    
    #     Fish$M<- log(0.01)/-as.numeric(Fish$MaxAge)
    
    Fish$M<- Fish$vbk*Fish$MvK
    
    if (is.na(Fish$MaxAge))
    {
      Fish$MaxAge<- -log(0.01)/Fish$M
    }
    
    #     if (is.na(SpeciesLifeHistory$Mat50)==F & is.na(SpeciesLifeHistory$LocalMaturity)==F)
    #     {
    #       Fish$Linf<- SpeciesLifeHistory$Mat50/0.7
    #     }
    #     
    
    if (is.na(SpeciesLifeHistory$Mat50) ) #Use Prince et al. 2014 LHI
    {
      Fish$Mat50<- Fish$Linf*Fish$LengthMatRatio
      
      Fish$Mat95<- as.numeric(1.05*Fish$Mat50)
    }
    
    if (is.na(SpeciesLifeHistory$AgeMat50)==F)
    {
      Fish$Mat50<- LengthAtAge(SpeciesLifeHistory$AgeMat50,Fish,0)
      
      Fish$Mat95<- as.numeric(1.05*Fish$Mat50)
    }
    
    Fish$Mat95<- as.numeric(1.05*Fish$Mat50)
    
    ReformData<- FormatCCFRPData(iGFD)
    
    LengthData<- ReformData$LengthData
    
    DensityData<- ReformData$DensityData
    
    write.csv(file=paste(Directory,AssessmentName,'_LengthData.csv',sep=''),LengthData)
    
    write.csv(file=paste(Directory,AssessmentName,'_DensityData.csv',sep=''),DensityData)
    
    FigureFolder<- paste(Directory,'Figures/',sep='')
    
    ResultFolder<- paste(Directory,'Results/',sep='')
    
    dir.create(FigureFolder,recursive=T)
    
    dir.create(ResultFolder,recursive=T)
    
    PlotLifeHistory()
    
    for (d in 1:length(AvailableData)) #Read in available data
    {
      #       eval(parse(text=paste(AvailableData[d],'<- read.csv(',"'",Directory,Fishery,'_',AvailableData[d],'.csv',"'",')',sep='')))
      eval(parse(text=paste('Plot',AvailableData[d],'(',AvailableData[d],')',sep='')))
    }
    
    ### Run Assessments ###
    
    AssessmentResults<- as.data.frame(matrix(NA,nrow=length(Assessments)*10,ncol=9))
    
    colnames(AssessmentResults)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
    
    AssessmentResults$Year<- as.numeric(AssessmentResults$Year)
    
    Count<-0
    
    Fish$LHITol<- 0.99
    
    # LengthData<- LengthData[LengthData$Year>2006,]
    
    for (a in 1:length(Assessments)) #Loop over possible assessments, store in Assessment results. Many assessments have more detailed outputs than can also be accessed 
    {
      
      if (Assessments[a]=='LBAR') #Run LBAR assessment
      {
        
        SampleCheck<- CheckLengthSampleSize(LengthData)        
        
        if (SampleCheck$YearsWithEnoughData>0)
        {
          
          #           LengthData<- SampleCheck$ParedData
          
          Temp<- LBAR(SampleCheck$ParedData,1,0.2,0,ReserveYear,NA,1000,1,1,NA)$Output		
          # Temp2<- OldLBAR(LengthData,1,0.2,0,100,1,1)$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp	
          
          Count<- Count+DataLength	
        }
      }
      
      if (Assessments[a]=='CatchCurve') #Run Catch Curve analysis
      {
        
        SampleCheck<- CheckLengthSampleSize(LengthData)        
        
        if (SampleCheck$YearsWithEnoughData>0)
        {
          
          Temp<- CatchCurve(SampleCheck$ParedData,'AgeBased',1,ReserveYear,NA,0,2,1,0,1)$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp	
          
          Count<- Count+DataLength
        }
      }
      
      
      if (Assessments[a]=='DensityRatio') #Run density ratio analysis 
      {
        Temp<- DensityRatio(DensityData,1,1,'Count',100,1)$Output
        
        ddply(DensityData,c('Year'),summarize,huh=length(Site))
        
        DataLength<- dim(Temp)[1]
        
        AssessmentResults[(Count+1):(Count+DataLength),]<- Temp	
        
        Count<- Count+DataLength
      }
      
      if (Assessments[a]=='CatchMSY')
      {
        Temp2<- CatchMSY(CatchData,3000,0.05,0,0,1,0,0,1,NA,c(0.75,0.99),NA,NA,c(0.25,0.65))
        
        Temp<- Temp2$Output
        
        DataLength<- dim(Temp)[1]
        
        AssessmentResults[(Count+1):(Count+DataLength),]<- Temp		
        
        Count<- Count+DataLength
        
      }
      
      if (Assessments[a]=='LBSPR') #Run LBSPR Assessment
      {
        
        SampleCheck<- CheckLengthSampleSize(LengthData)        
        
        if (SampleCheck$YearsWithEnoughData>0)
        {
          
          #           LengthData<- SampleCheck$ParedData
          #           SampleCheck$ParedData$Year[  SampleCheck$ParedData$Year<=2010]<- 1
          #           
          #           SampleCheck$ParedData$Year[  SampleCheck$ParedData$Year>=2011]<- 2
          
          Temp2<- LBSPR(SampleCheck$ParedData,0,10,1,1,1,ReserveYear)
          
          Temp<- Temp2$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp		
          
          Count<- Count+DataLength
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
        
        DataLength<- dim(Temp)[1]
        
        AssessmentResults[(Count+1):(Count+DataLength),]<- Temp		
        
        Count<- Count+DataLength
      }
      
      
      show(paste('Finished ',Assessments[a],'-',round(100*a/length(Assessments),2),'% Done',sep=''))
      
    }
    
    AssessmentResults<- AssessmentResults[is.na(AssessmentResults$Year)==F,]
    AssessmentResults$Year<- as.numeric(AssessmentResults$Year)
    AssessmentResults$Value<- as.numeric(AssessmentResults$Value)
    AssessmentResults$LowerCI<- as.numeric(AssessmentResults$LowerCI)
    AssessmentResults$UpperCI<- as.numeric(AssessmentResults$UpperCI)
    AssessmentResults$SD<- as.numeric(AssessmentResults$SD)
    
    AssessmentResults[,4:7]<- round(AssessmentResults[,4:7],2)
    
    show(AssessmentResults)
    
    if (sum(AssessmentResults$Method=='CatchCurve')>0)
    {
      
      SummaryPanel(AssessmentResults,LengthData,Species,Sites[s],YearsToSmooth=3)
      
    }
    save.image(file=paste(ResultFolder,AssessmentName,'_Settings.RData',sep='')) #Save settings used to produce current results
    write.csv(file=paste(ResultFolder,AssessmentName,'_Results.csv',sep=''),AssessmentResults) #Save current results
    
  } #Close species  (f) loop
} #Close sites (s) loop




