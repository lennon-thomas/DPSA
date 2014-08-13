FormatCCFRPData<- function(Data)
{
  # Format Length Data ------------------------------------------------------
  
  #   Data<- iGFD
  
  LengthDataNames<- c('Year','Month','Site','Length','LengthType','Sex','Special','MPA','FisheryDependent')
  
  LengthData<- as.data.frame(matrix(NA,nrow=dim(Data)[1],ncol=length(LengthDataNames)))
  
  colnames(LengthData)<- LengthDataNames
  
  LengthData$Year<- Data$Year
  
  LengthData$Month<- Data$Month
  
  LengthData$Site<- Data$Site
  
  LengthData$Length<- Data$length_cm
  
  LengthData$LengthType<- 'cm'
  
  LengthData$Sex<- 'Unknown'
  
  LengthData$Special<- paste('Gear is ',Data$Sample_Type,sep='')
  
  LengthData$MPA<- Data$MPA_or_REF
  
  LengthData$MPA[LengthData$MPA=='REF']<- 0
  
  LengthData$MPA[LengthData$MPA=='MPA']<- 1
  
  LengthData$MPA<- as.numeric( LengthData$MPA)
  
  LengthData$FisheryDependent<- 1
  
  # Format Density Data -----------------------------------------------------
  
  Data$Weight<- Fish$WeightA* Data$length_cm ^ Fish$WeightB
    
  DensityData<- ddply(Data,c('Year','Month','sample_Idcellday'),summarize,Site='All',Count=length(length_cm),Biomass=sum(Weight,na.rm=T)
                      ,SampleArea= mean(Sample_Area,na.rm=T),AreaUnits=unique(Area_units),DistanceFromBorder=mean(Meters.to.MPA.border,na.rm=T)
                      ,SampleType=unique(Sample_Type),MPA=unique(MPA_or_REF))
  
  DensityData$MPA[DensityData$MPA=='REF']<- 0
  
  DensityData$MPA[DensityData$MPA=='MPA']<- 1
  
  return(list(LengthData=LengthData,DensityData=DensityData))
  
}




