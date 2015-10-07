#### ControlFile ###

AvailableData<- c("LengthData")

# Assessments<- c('CatchCurve')

 Assessments<- c('LBAR','LBSPR')
 
 Fish<- NULL
 
 Fish$SciName<- 'O. '
 
 Fish$CommName<- 'Yellowtail snapper'
### Life History ###
 
 DefaultSD<- 0.05
 
 Fish$LHITol<- 0.1 # The average % deviation from LHI allowed
 
 Fish$vbk<- 0.133                      ##DO:according to Informe final   0.21625    error 0.004M and 0.06H
 
 Fish$LengthError<- 0.05
 
 Fish$Linf<- 49.974                       ##41.9M and 40.14H averaged from Norahs values (used averages = 41.02) used our equation to determine Linf for tail length (average = 25.15 and H = 24.57)
 
 Fish$t0<- -3.132
 
 Fish$WeightA <-0.000117
 
 Fish$WeightB<- 2.6504
 
 

 
 
 Fish$M<- 0.194
 
 
 Fish$MvK<- 1.46 #Ratio of M versus K
 
 #Fish$MaxAge<-  #Max age is age at thich only 1% of population are left. Can also be set manually if needed
 
 Fish$MortalityError<- 0.1                 ## test error, can change to accomodate more error
 
 Fish$Mat50<- 27.758                      ## 24.05 according to Toral and Hearn (tail only = 13.96)
 
 Fish$Mat95<- 31.644                              ## 26 according to Toral and Hearn (tail only = 15.24)
 
 Fish$PLD<- NA
 
 Fish$VBSD<- 0.05
 
 Fish$VBErrorSlope<- 0.1
 
 Fish$res<- 'Low'
 
 ### Fleet Parameters ###
 
 Fleet<- NULL
 
 Fleet$MinSizeCaught<- 17                        ##26 cm minimum landing size for Galapagos      (was 20 before?) (tail = 15.24) 
 
 Fleet$MaxSizeCaught<- Fish$Linf-0.5
 
