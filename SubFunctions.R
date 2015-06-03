
PlotLifeHistory<- function()
{
  
  Lengths<- LengthAtAge(1:Fish$MaxAge,Fish,0)
  
  Maturity<- MaturityAtAge(Lengths,Fish)
  
  LifeMat<- data.frame(1:Fish$MaxAge,Lengths,Maturity)
  
  colnames(LifeMat)<- c('Age','Length','PercentMature')
  
  pdf(file=paste(FigureFolder,AssessmentName,' Life History.pdf',sep=''))
  
  LifePlot<- ggplot(LifeMat,aes(x=Age,y=Length,color=PercentMature))+geom_line(size=2)+ylab('Length (cm)')
  
  print( LifePlot + scale_colour_gradient(limits=c(0, 1), low="steelblue2",high='red'))  
  
  dev.off()
}



PlotCatchData<- function(CatchDat)
{
  
  library(zoo)
  # CatchDat=CatchData
  
  InterpCatch<- na.approx(CatchDat$Catch)
  
  InterpMarker=as.numeric(is.na(CatchDat$Catch))
  
  PointStyle=21*InterpMarker
  PointStyle[PointStyle==0]<- 16
  
  pdf(file=paste(FigureFolder,' Catch History.pdf',sep=''),width=7,height=4) 
  par(mai=c(1,1,1,1.5))
  Month<- CatchDat$Month
  Month[Month==-999]<- 'Total'
  TimeName=paste(CatchDat$Year,Month,sep='-')
  plot(InterpCatch,xaxt='n',pch=PointStyle,col=InterpMarker+1,ylab=CatchDat$Units[1],xlab=NA,cex=1.2,pty='m',bty='n')
  axis(1,at=1:length(TimeName),labels=TimeName,las=0)
  legend('right',pch=c(16,21),col=c(1,2),legend=c('Real','Interpolated'),xpd=T,bty='n',inset=-.3)
  dev.off()
  
}


ApplyLifeHistoryError<- function()
{
  
#   LHI_Error <- 1.1
  
#   while(LHI_Error>Fish$LHITol) # Only accept values that don't violate LHI too much
#   {
    NewFish<- Fish
    
    NewFish$vbk<- rlnorm(1,log(Fish$vbk),Fish$LengthError)
    
    NewFish$MvK<- runif(1,NewFish$MinMvK,NewFish$MaxMvK)
    
    NewFish$M<- NewFish$vbk * NewFish$MvK
    
#     NewFish$M<-rlnorm(1,log(Fish$M),(Fish$MortalityError))
    
    NewFish$Linf<- rlnorm(1,log(Fish$Linf),Fish$LengthError)
    
    NewFish$LengthMatRatio<- runif(1,NewFish$MinLengthMatRatio,NewFish$MaxLengthMatRatio)

    NewFish$Mat50<- NewFish$Linf * NewFish$LengthMatRatio

    NewFish$Mat95<- NewFish$Mat50*1.05    

#     if (Fish$t0!=0)
#     {
#       NewFish$t0<- rnorm(1,(Fish$t0),(Fish$LengthError))
#     }
#       
    
    NewFish$MaxAge<- -log(.01)/Fish$M
        
#     Mtm_Error<- abs((NewFish$M*AgeAtLength(NewFish$Mat95,Fish,0))/1.65-1)
#     
#     MvK_Error<- abs((NewFish$M/NewFish$vbk)/1.6-1)
#     
#     LmvLinf_Error<- abs((mean(c(NewFish$Mat50,NewFish$Mat95))/NewFish$Linf)/0.67-1)
#     
#     LHI_Error<- (mean(c(Mtm_Error,MvK_Error, LmvLinf_Error)))
    
#   }
  
  return(NewFish)
}

LBSPR_SingleSpeciesAssessmentfun<- function(CatchatLength,AssessDir,CurrentDir,LengthBins,Year,EstimatedM,Fish)
{
  Output<- as.data.frame(matrix(NA,nrow=1,ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  Flag<- 'None'
  # Year=Years[y]
  Details<- as.data.frame(matrix(NA,nrow=1,ncol=4))
  
  colnames(Details)<- c('Year','FvM','SelL50','SelL95')
  #Modified from the version sent by Sarah V. in Sep 2013
  ##############################
  # Read in Assumed Parameters #
  ##############################
  
  SpeciesName <- Species
  # assumedLinf  <- 1.2*max(CatchatLength) 
  assumedLinf  <- Fish$Linf
  
  M<- Fish$M
  
  if (is.numeric(EstimatedM)==T)
  {
    M<- EstimatedM
  }	
  
  assumedMK  <- Fish$MvK
  
  # assumedMK  <- M/Fish$vbk
  
  
  genM    <-   0.1 # as.numeric(as.character(SpFile[4, 1 +1]))
  genLinf		<- 1
  genLinfCV	<- Fish$LengthError
  gent0		<- Fish$t0
  MatType		<- 'Logistic'
  MatL50		<- Fish$Mat50  
  MatL95		<- Fish$Mat95
  Wbeta		<- Fish$WeightB
  Walpha		<- Fish$WeightA
  genK    <- genM/assumedMK
  # minLen  <- Fleet$MinSizeCaught
  minLen  <- min(CatchatLength,na.rm=T)
  maxLen  <- 1.2 * max(CatchatLength,na.rm=T)
  # maxLen  <- (Fleet$MaxSizeCaught)
  # maxLen  <- 1.1*max(CatchatLength)    
  #################################################################
  # Convert Data received to appropriate format        
  #################################################################
  
  LengthClasses <- seq(floor(minLen),ceiling(maxLen),by = LengthBins)
  LengthMids <- seq(LengthClasses[1] +((LengthClasses[2]-LengthClasses[1])/2), by=(LengthClasses[2]-LengthClasses[1]), length=length(LengthClasses)-1)
  # LenFreq <- hist(CatchatLength,breaks=seq(floor(minLen),ceiling(maxLen),by=LengthBins),plot=FALSE,right=F)$counts
  LenFreq<- DanHist(CatchatLength,seq(floor(minLen),ceiling(maxLen),by=LengthBins))$Frequency
  
  # DanHist(CatchatLength,seq(minLen,maxLen+LengthBins,by=LengthBins))
  
  LenProp <- as.vector(LenFreq/sum(LenFreq))
  
  ########################################
  # Set final params to pass to ADMB file#    
  ########################################
  
  MK <- assumedMK
  LinfTest <- assumedLinf
  LinfCV <- genLinfCV
  PercLeft <- 0.01
  NumAgClass <- 100
  L50 <-  MatL50
  L95 <-  MatL95
  
  ###########################################
  # Set working directory and run assessmemt#    
  ###########################################
  
  #MNew Assess Dir for each method tested
  setwd(AssessDir)
  WriteDat(MK, LinfTest, LinfCV, PercLeft, NumAgClass, LengthMids, LengthClasses, LenFreq, LenProp, L50, L95, AssessDir) 
  
  #Find best Starting Values
  LengthComp <- cbind(LengthMids, LenFreq)
  # Write Pin
  FirstLen <- LengthComp[min(which(LengthComp[,2] > 0)),1]
  # Find best starting values
  Vals1  <- c(FirstLen/LinfTest, 0.05, log(1))
  Vals2 <- c(FirstLen/LinfTest+0.1, 0.05, log(0.5))
  Vals3 <- c(FirstLen/LinfTest+0.1, 0.05, log(2.5))
  
  ValMat <- matrix(c(Vals1, Vals2, Vals3), byrow=T, nrow=3)
  
  tempRunAssess <-function (Vals) {
    WritePin(AssessDir, Vals)
    ADMBFile <- paste(AssessDir, "/LBSPR_AssessFun", sep="")
    ModelFailed <- FALSE
    setwd(AssessDir)
    #system(ADMBFile)
    run_admb("LBSPR_AssessFun") 
    # Check Model Failed 
    FileList <- c("admodel.cov", "lbspr_assessfun.cor", "lbspr_assessfun.std")
    allFiles <- list.files(AssessDir)
    if(any(FileList %in% allFiles) == FALSE) ModelFailed <- TRUE
    if (ModelFailed) {
      estSel50 <- NA
      estSel95 <- NA
      estFM <- NA
      ObjVal <- 1E6
      estSPR <- NA
      ModelFit <- NA
    } else {
      estSel50 <- read.table(paste(AssessDir, "/LBSPR_AssessFun.par", sep=""))[1,1] * LinfTest
      estSel95 <- estSel50 + read.table(paste(AssessDir, "/LBSPR_AssessFun.par", sep=""))[2,1] * LinfTest
      estFM <- exp(read.table(paste(AssessDir, "/LBSPR_AssessFun.par", sep=""))[3,1])
      ObjVal <- scan(paste(AssessDir, "/LBSPR_AssessFun.rep", sep=""), what=double(), skip=2, nlines=1, quiet=TRUE)
      estSPR <- scan(paste(AssessDir, "/LBSPR_AssessFun.rep", sep=""), what=double(), skip=3, nlines=1, quiet=TRUE)
      ModelFit <- scan(paste(AssessDir, "/LBSPR_AssessFun.rep", sep=""), what=double(), nlines=1, quiet=TRUE) 
    }
    setwd(CurrentDir)
    Output <- NULL
    Output$est <- c(estSel50, estSel95, estFM, estSPR, ObjVal) 
    Output$fit <- ModelFit
    return(Output)
  }  
  
  DeleteFiles(AssessDir) 
  SaveResults <- matrix(NA, nrow=nrow(ValMat), ncol=5)
  SaveFit <- rep(list(NA), nrow(ValMat))
  for (i in 1:nrow(ValMat)) {
    temp <- tempRunAssess(ValMat[i,]) 
    SaveResults[i,] <- temp$est
    SaveFit[[i]] <- temp$fit
  }
  
  ModelFailed <- FALSE
  if (all(is.na(SaveResults[,1:4]))) ModelFailed <- TRUE
  
  MinInd <- which.min(SaveResults[,5])
  
  ###########################################
  # Unpack parameters#    
  ###########################################
  
  LBSPR_Output <- NULL
  LBSPR_Output$SelL50 <- SaveResults[MinInd, 1]
  LBSPR_Output$SelL95 <- SaveResults[MinInd, 2]
  LBSPR_Output$EstFM <- SaveResults[MinInd, 3]
  LBSPR_Output$EstSPR <- SaveResults[MinInd, 4]
  LBSPR_Output$ObjVal <- SaveResults[MinInd, 5]
  LBSPR_Output$ModelFit <- SaveFit[[MinInd]]
  LBSPR_Output$ModFAILED <- ModelFailed
  
  Output<-NULL
  
  Output$Year<- Year
  
  Output$Method<- 'LBSPR'
  
  Output$Value<- LBSPR_Output$EstSPR
  
  Output$Metric<- 'SPR'
  
  Output$LowerCI<- NA
  
  Output$UpperCI<- NA
  
  Output$SD<- NA
  
  Output$Flag<- paste('ModelFailed is',ModelFailed)
  
  Details <- NULL
  
  Details$Year<- Year
  
  Details$FvM <- LBSPR_Output$EstFM
  
  Details$SelL50<- LBSPR_Output$SelL50
  
  Details$SelL95<- LBSPR_Output$SelL95
  
  # Do Plot function
  DoSinglePlot <- function(ModelFit, LenProb, LenMids, SaveFile=TRUE, FileName, ModFail) {
    
    pdf(file=paste(CurrentDir,'/',FigureFolder,FileName,sep=''))
    par(mfrow=c(1,1))
    Max <- max(LenProb)
    xx <- barplot(LenProb, names.arg=round(LenMids,2), ylim=c(0, Max+Max*0.1))
    if (ModFail == FALSE) lines(xx, ModelFit, lwd=4)
    if (ModFail) text(xx[5], Max*0.8, "MODEL FAILED TO CONVERGE", cex=1.5, pos=4)
    mtext(side=1, "Length Classes", line=3, cex=1.5)
    mtext(side=2, "Proportion", line=2.5, cex=1.5)
    dev.off()
  }
  
  DoSinglePlot(LBSPR_Output$ModelFit,LenProp,LengthMids,SaveFile=TRUE,FileName=paste(Year,' LBSPRModelFit.pdf'),ModelFailed)
  
  
  
  LengthMids[LengthMids>Fish$Linf]<- 0.99*Fish$Linf
  AgeVector<- floor(AgeAtLength(LengthMids,Fish,0))
  
  Ages<- unique(AgeVector[is.na(AgeVector)==F])
  
  Residuals<- LenProp - LBSPR_Output$ModelFit
  
  
  CohortDeviates<- as.data.frame(matrix(NA,nrow=length(Ages),ncol=2))
  
  AgeDeviates<- as.data.frame(matrix(NA,nrow=length(Ages),ncol=2))
  
  
  for (a in 1:length(Ages))
  {
    Where<- (AgeVector==Ages[a] & is.na(AgeVector)==F)
    
    CohortDeviates[a,]<- data.frame(Year-Ages[a],sum(Residuals[Where]))
    
    AgeDeviates[a,]<- data.frame(Ages[a],sum(Residuals[Where]))
    
  }
  
  colnames(CohortDeviates)<- c('Cohort','Residuals')
  
  colnames(AgeDeviates)<- c('Age','Residuals')
  
  
  
  if (ModelFailed== TRUE)
  {
    Output$Year<- Year
    
    Output$Method<- 'LBSPR'
    
    Output$Flag<- paste('ModelFailed is',ModelFailed)
    
  }
  
  return(list(Output=Output,Details=Details,CatchCurveResiduals= CohortDeviates,AgeResiduals= AgeDeviates))
}

LengthAtAge<- function(Ages,Fish,Error)
{
  
  LenSD<- Error*(1+Fish$VBErrorSlope*Ages/Fish$MaxAge)
  
  RawLengths<- Fish$Linf*(1-exp(-Fish$vbk*(Ages-Fish$t0)))
  
  LengthWithError<- RawLengths*rlnorm(length(Ages),mean=0,sd=LenSD)
  
  return(LengthWithError)
}

AgeAtLength<- function(Lengths,Fish,Error)
{
  # Error<- Fish$LengthError	
  Lengths[is.na(Lengths)]<- 0
  # Lengths<- LengthDat$Length
  AgeSD<- Error*(1+Fish$VBErrorSlope*Lengths/Fish$Linf)
  #   RawAges<- (log(1-(Lengths)/Fish$Linf)/-Fish$vbk)+Fish$t0
  RawAges<- (log(1-pmin(Lengths,Fish$Linf*.99)/Fish$Linf)/-Fish$vbk)+Fish$t0
  AgeWithError<- RawAges*rlnorm(length(Lengths),mean=0,sd=AgeSD)
  
  return(AgeWithError)
}

DanHist<- function(Data,Breaks)
{
  # Data<- TempLengthDat$Age[TempLengthDat$MPA==1]
  
  # Breaks<- BinBreaks
  
  BreakStore<- as.data.frame(matrix(NA,nrow=length(Breaks),ncol=3))
  colnames(BreakStore)<- c('Age','Frequency','LogFrequency')
  for (b in 1:(length(Breaks)-1))
  {
    BreakStore[b,1:2]<- c(Breaks[b],sum(Data>= Breaks[b] & Data< Breaks[b+1],na.rm=T))
    BreakStore[b,3]<- log(BreakStore[b,2])
    
  }
  
  BreakStore$LogFrequency[is.infinite(BreakStore$LogFrequency)]<- NA
  return(BreakStore[1:(length(Breaks)-1),])
}


CalculateDensity<- function(Densities,Years,Weights,Form)
{
  
  #   Years<- LaggedYears
  #   Densities<- TempDenDat
  #   
  #   # # # Densities<- DenDat[DenDat$Year %in% Years]	
  #   Form<- 'Biomass'
  #   Weights<- weights
  
  Densities$DistanceFromBorder[Densities$DistanceFromBorder==-999]<- NA
  
  Densities$DistanceFromBorder[is.na(Densities$DistanceFromBorder)]<- mean(Densities$DistanceFromBorder,na.rm=T)
  
  
  DensityForm<- colnames(Densities)==Form
  
  LagDensity<- as.data.frame(matrix(NA,nrow= length(Years),ncol=4))
  
  colnames(LagDensity)<- c('Year','MPADensity','FishedDensity','DensityRatio')
  
  WeightedDensity<- as.data.frame(matrix(NA,nrow=1,ncol=4))
  
  colnames(WeightedDensity)<- c('Year','MPADensity','FishedDensity','DensityRatio')	
  for (y in 1:length(Years))
  {
    
    YearlyDensity<- Densities[Densities$Year==Years[y],]
    
    Reserve<- YearlyDensity$MPA==1
    
    MPADensity<- sum(YearlyDensity$DistanceFromBorder[Reserve]*(YearlyDensity[Reserve,DensityForm]/YearlyDensity$SampleArea[Reserve]))/sum(YearlyDensity$DistanceFromBorder[Reserve])
    
    FishedDensity<- sum(YearlyDensity$DistanceFromBorder[Reserve==F]*(YearlyDensity[Reserve==F,DensityForm]/YearlyDensity$SampleArea[Reserve==F]),na.rm=T)/sum(YearlyDensity$DistanceFromBorder[Reserve==F],na.rm=T)
    
    #     MPADensity<- sum(YearlyDensity$DistanceFromBorder[Reserve]*YearlyDensity[Reserve,DensityForm])/sum(YearlyDensity$DistanceFromBorder[Reserve]*YearlyDensity$SampleArea[Reserve])
    #     
    #     FishedDensity<- sum(YearlyDensity$DistanceFromBorder[Reserve==F]*YearlyDensity[Reserve==F,DensityForm],na.rm=T)/sum(YearlyDensity$DistanceFromBorder[Reserve==F]*YearlyDensity$SampleArea[Reserve==F],na.rm=T)
    #     
    LagDensity[y,]<- data.frame(Years[y],MPADensity,FishedDensity,FishedDensity/MPADensity)
    
  }
  WeightedDensity[1,]<- data.frame(Years[length(Years)],sum(Weights*LagDensity$MPADensity)/sum(Weights),sum(Weights*LagDensity$FishedDensity)/sum(Weights),sum(Weights*LagDensity$DensityRatio)/sum(Weights))
  
  return(WeightedDensity)
  
}

CalculateCPUE<- function(CPUE,Years,Weights,Form)
{
  
  CPUE$DistanceFromBorder[CPUE$DistanceFromBorder==-999]<- NA
  
  CPUE$DistanceFromBorder[is.na(CPUE$DistanceFromBorder)]<- mean(CPUE$DistanceFromBorder,na.rm=T)
  
  
  CPUEForm<- colnames(CPUE)==Form
  
  LagCPUE<- as.data.frame(matrix(NA,nrow= length(Years),ncol=4))
  
  colnames(LagCPUE)<- c('Year','MPACPUE','FishedCPUE','CPUERatio')
  
  WeightedCPUE<- as.data.frame(matrix(NA,nrow=1,ncol=4))
  
  colnames(WeightedCPUE)<- c('Year','MPACPUE','FishedCPUE','CPUERatio')
  for (y in 1:length(Years))
  {
    
    YearlyCPUE<- CPUE[CPUE$Year==Years[y],]
    
    Reserve<- YearlyCPUE$MPA==1
    
    MPACPUE<- sum(YearlyCPUE$DistanceFromBorder[Reserve]*(YearlyCPUE[Reserve,CPUEForm]/YearlyCPUE$AnglerHours[Reserve]))/sum(YearlyCPUE$DistanceFromBorder[Reserve])
    
    FishedCPUE<- sum(YearlyCPUE$DistanceFromBorder[Reserve==F]*(YearlyCPUE[Reserve==F,CPUEForm]/YearlyCPUE$AnglerHours[Reserve==F]),na.rm=T)/sum(YearlyCPUE$DistanceFromBorder[Reserve==F],na.rm=T)
    
    LagCPUE[y,]<- data.frame(Years[y],MPACPUE,FishedCPUE,FishedCPUE/MPACPUE)
    
  }
  WeightedCPUE[1,]<- data.frame(Years[length(Years)],sum(Weights*LagCPUE$MPACPUE)/sum(Weights),sum(Weights*LagCPUE$FishedCPUE)/sum(Weights),sum(Weights*LagCPUE$CPUERatio)/sum(Weights))
  
  return(WeightedCPUE)
  
}


movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  return(s/count)
}

MaturityAtAge <- function(Length,Fish)
{
  
  
  s50<- Fish$Mat50
  
  s95<- Fish$Mat95
  
  mature<- ((1/(1+exp(-log(19)*((Length-s50)/(s95-s50))))))  
}


