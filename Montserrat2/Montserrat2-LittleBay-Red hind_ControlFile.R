#### ControlFile ###

AvailableData<- c("LengthData")

# Assessments<- c('CatchCurve')

 Assessments<- c('LBSPR')
 
 Fish<- NULL
 
 Fish$SciName<- 'E. guttas'
 
 Fish$CommName<- 'Red hind'

### Life History ###
 
 DefaultSD<- 0.05
 
 Fish$LHITol<- 0.1 # The average % deviation from LHI allowed
 
 Fish$vbk<- 0.0859                       ##DO:according to Informe final   0.21625    error 0.004M and 0.06H
 
 Fish$LengthError<- 0.05
 
 Fish$Linf<- 55.78                        ##41.9M and 40.14H averaged from Norahs values (used averages = 41.02) used our equation to determine Linf for tail length (average = 25.15 and H = 24.57)
 
 Fish$t0<- -3.817
 
 Fish$WeightA <- 0.0134 #From D. Viana's summary of lit, "Lobster Life History..."
 
 Fish$WeightB<- 3.0402#From D. Viana's summary of lit, "Lobster Life History..."
 
 
 # Fish$MaxAge<- 25
 
 #Fish$M<- 0.17                   ## average mortality according to Matt Kay and literature review, Hearn 2008 use 0.348M and 0.378H (average = 0.363)
 
 Fish$M<-    .171185
 
 # Fish$M<- 0.166  
 
 
 Fish$MvK<- 1.992841 #Ratio of M versus K
 
 Fish$MaxAge<- log(.01)/-Fish$M #Max age is age at thich only 1% of population are left. Can also be set manually if needed
 
 Fish$MortalityError<- 0.1                 ## test error, can change to accomodate more error
 
 Fish$Mat50<- 24.70                      ## 24.05 according to Toral and Hearn (tail only = 13.96)
 
 Fish$Mat95<- 28.18                              ## 26 according to Toral and Hearn (tail only = 15.24)
 
 Fish$PLD<- NA
 
 Fish$VBSD<- 0.05
 
 Fish$VBErrorSlope<- 0.1
 
 Fish$res<- 'Low'
 
 ### Fleet Parameters ###
 
 Fleet<- NULL
 
 Fleet$MinSizeCaught<- 22                        ##26 cm minimum landing size for Galapagos      (was 20 before?) (tail = 15.24) 
 
 Fleet$MaxSizeCaught<- Fish$Linf-0.5
 
