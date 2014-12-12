
# AddMissingFish ----------------------------------------------------------
#This code adds "zero" counts for all species ever observed at a site but not observed on a given trip
AddMissingFish<- function(Data)
{
  
#   Data<- GFD
  
  SpeciesTable<- unique(Data [,colnames(Data)[16:dim(Data)[2]]])
  
  SpeciesSightings<- ddply(Data,c('Site'),summarize,SpeciesSeen=unique(CommName))
  
  SpeciesSightingsByTrip<- ddply(Data,c('sample_Idcellday'),summarize,SpeciesSeen=length(unique(CommName)))
  
  BlankVars<- colnames(Data)
  
  BlankVars<- BlankVars[!(BlankVars%in%c('Year','Month','Site','sample_Idcellday','Sample_Type','Sample_Area','Area_units',
                                         'Angler_hours','MPA_or_REF','GRID_ID_cell','Meters.to.MPA.border'))]
  Sites<- unique(Data$Site)
  
  for (s in 1:length(Sites))
  {
    
    Trips<- unique(Data$sample_Idcellday[Data$Site==Sites[s]])
    
    SpeciesSeen<- SpeciesSightings$SpeciesSeen[SpeciesSightings$Site==Sites[s]]
    
    for (t in 1:length(Trips))
    {
      
      SpeciesSpotted<- unique(Data$CommName[Data$sample_Idcellday==Trips[t]])
      
      SpeciesMissing<- SpeciesSeen[!(SpeciesSeen %in% SpeciesSpotted)]
      
      if (sum(!(SpeciesSeen %in% SpeciesSpotted))>0)
      {
        BlankTrip<- Data[Data$sample_Idcellday==Trips[t],][1,]
        
        BlankTrip<- RepMat(BlankTrip,length(SpeciesMissing),'Rows')
        
        BlankTrip[,colnames(BlankTrip) %in% BlankVars]<- NA
        
        MissingData<- SpeciesTable[SpeciesTable$CommName %in% SpeciesMissing,]
        
        BlankTrip[,'length_cm']<- 0
        
        BlankTrip [,colnames(Data)[16:dim(Data)[2]]]<- MissingData
        
        Data<- rbind(Data,BlankTrip)
        
      } #Close if any missing loop
      
    } #Close Trips
    
  } #Close sites
  
  Data<- Data[order(Data$Year,Data$Month,Data$Site,Data$sample_Idcellday),]
  
  return(Data)
} #Close function
