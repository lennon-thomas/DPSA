###### CCFRP Data Poor Stock Assessments ######
#This wrapper calls the AssessmentModules and SubFunctions of the DPSA modules 
#created by Ovando et al. to run DPSAs on the CCFRP groundfish data

#Created by Dan Ovando


### Setup Working Environment ###
rm(list=ls())
set.seed(442)

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
library(grid)

sapply(list.files(pattern="[.]R$", path="Functions", full.names=TRUE), source)

# source("AssessmentModules.R") #Pull in assessment modules
source("SubFunctions.R") #Pull in helper functions for assessment modulesAssessmentName<- 'LBSPR'

Country<- 'Montserrat'

Site<- 'LittleBay'

Species<- 'Red hind'

Fishery <- paste(Country,Site,Species,sep='-')

Directory<- paste(Country, "/", Site,'/',Species,'/',sep='')

Font <- 'Helvetica'

FontColor <- 'Black'

PlotFontSize <- 11
MPAColor <- "#1b9e77"

FishedColor <- "#d95f02"
Theme<- theme(legend.position='top',plot.background=element_rect(color=NA),
              rect=element_rect(fill='transparent',color=NA)
              ,text=element_text(size=12,family=Font,color=FontColor),
              axis.text.x=element_text(color=FontColor),
              axis.text.y=element_text(color=FontColor),legend.key.size=unit(1,'cm'))


source(paste(Country, "/", Site,'/',Species,'/',Fishery,"_ControlFile.R", sep = ""))


dir.create(paste(Directory,'Results',sep=''))

dir.create(paste(Directory,'Figures',sep=''))

FigureFolder<- paste(Directory,'Figures/',sep='')

ResultFolder<- paste(Directory,'Results/',sep='')

for (d in 1:length(AvailableData)) #Read in available data
{
  eval(parse(text=paste(AvailableData[d],'<- read.csv(',"'",Directory,Fishery,'_',AvailableData[d],'.csv',"'",')',sep='')))
  
  eval(parse(text=paste('Plot',AvailableData[d],'(LengthData,FigureFolder,Fish,Species,Site,Theme)',sep='')))
  
}
Fish$AgeMat50<- NA

Fish$AgeMatSource<- NA

AssessmentResults<- list()

MonteResults<- list()

Counter<- 0

Fishes<- Species
#############################
###########LBAR##############
############################

Temp<- LBAR(LengthData,LagLength=1,Weight=1,IncludeMPA=0,ReserveYr=NA,OutsideBoundYr=NA,Iterations=100,
            BootStrap=1,LifeError=1,Lc=NA)$Output		

#############################################################################################################
## Run LBSPR
#############################################################################################################
test<-"C:/Users/lthomas/GitHub/DPSA/"
MinSampleSize<-15
NumIterations<-1

SampleCheck<- CheckLengthSampleSize(LengthData)        

if (SampleCheck$YearsWithEnoughData>0)
{
  
  LengthQuantile<- quantile(SampleCheck$ParedData$Length,na.rm=T)
  
  
  #           Temp2<- LBSPR(SampleCheck$ParedData,EstimateM=0,Iterations=1,BootStrap=1,
  #                         LifeError=1,LengthBins=1,ReserveYear=ReserveYear,SL50Min=LengthQuantile[1],
  #                         SL50Max=LengthQuantile[2],DeltaMin=NA,DeltaMax=NA,IncludeReserve=TRUE)
 ## CatchatLength,AssessDir,CurrentDir,LengthBins,Year,EstimatedM,Fish
  Temp2<- LBSPR_SingleSpeciesAssessmentfun(SampleCheck$ParedData, Directory,test, LengthBins=1,Year=2015,EstimatedM,Fish)
  
  MonteCarlo<- Temp2$Details
  
  StoreMonte<- data.frame(Species,Sites[s],Assessments[a],MonteCarlo[,c('Iteration','Year','FvM','SPR')],stringsAsFactors=F) %>%
    gather('Metric','Value',FvM,SPR,convert=T) %>%
    rename(Site=Sites.s.,Assessment=Assessments.a.) %>%
    select(Species,Site,Assessment,Iteration,Year,Value,Metric)
  
  StoreMonte$Metric[StoreMonte$Metric=='FvM']<- 'F/M (LBSR)'
  
  MonteResults[[Counter]]<- StoreMonte
  
  Temp<- Temp2$Output
  
  StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
    rename(Site=Sites.s.,Assessment=Assessments.a.)
  
  AssessmentResults[[Counter]]<-StoreAssess
}
###############################
########'LBSPR'################
##############################
MinSampleSize<-15


if (Assessments[a]=='LBSPR') #Run LBSPR Assessment
{
  
  SampleCheck<- CheckLengthSampleSize(LengthData)        
  
  if (SampleCheck$YearsWithEnoughData>0)
  {
    
    LengthQuantile<- quantile(SampleCheck$ParedData$Length,na.rm=T)
    
    
    #           Temp2<- LBSPR(SampleCheck$ParedData,EstimateM=0,Iterations=1,BootStrap=1,
    #                         LifeError=1,LengthBins=1,ReserveYear=ReserveYear,SL50Min=LengthQuantile[1],
    #                         SL50Max=LengthQuantile[2],DeltaMin=NA,DeltaMax=NA,IncludeReserve=TRUE)
   ( CatchatLength,AssessDir,CurrentDir,LengthBins,Year,EstimatedM,Fish)
    Temp2<-LBSR_SingleSpecies(SampleCheck$ParedData,EstimateM=0,Iterations=1,BootStrap=1,
                  LifeError=1,LengthBins=1,ReserveYear=NA,SL50Min=LengthQuantile[1],
                  SL50Max=LengthQuantile[2],DeltaMin=0.01,DeltaMax=.5*Fish$Linf,IncludeReserve=FALSE)

    MonteCarlo<- Temp2$Details
    
    StoreMonte<- data.frame(Species,Sites[s],Assessments[a],MonteCarlo[,c('Iteration','Year','FvM','SPR')],stringsAsFactors=F) %>%
      gather('Metric','Value',FvM,SPR,convert=T) %>%
      rename(Site=Sites.s.,Assessment=Assessments.a.) %>%
      select(Species,Site,Assessment,Iteration,Year,Value,Metric)
    
    StoreMonte$Metric[StoreMonte$Metric=='FvM']<- 'F/M (LBSR)'
    
    MonteResults[[Counter]]<- StoreMonte
    
    Temp<- Temp2$Output
    
    StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
      rename(Site=Sites.s.,Assessment=Assessments.a.)
    
    AssessmentResults[[Counter]]<-StoreAssess
  }
}

