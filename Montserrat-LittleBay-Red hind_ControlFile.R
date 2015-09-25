#### ControlFile ###

AvailableData<- c("LengthData")

# Assessments<- c('CatchCurve')

 Assessments<- c('LBSPR')

### Life History ###

Fish<- NULL

Fish$vbk<- 0.0859

DefaultSD<- 0.05

Fish$Linf<- 55.78 #ish, 120mm SL

Fish$LengthError<- DefaultSD

Fish$t0<- -3.817

Fish$WeightA <- 0.0111
Fish$WeightB<- 3.1124

Fish$M<- .171185

Fish$MaxAge<- log(.01)/-Fish$M #Max age is age at thich only 1% of population are left. Can also be set manually if needed

Fish$MortalityError<- DefaultSD

Fish$Mat50<- 21.50

Fish$Mat95<- 24.51

#Fish$PLD<- 31
Fish$MvK<- NA


Fish$VBSD<- 0.03

Fish$VBErrorSlope<- 0.1
FecB   <- 3.1124


### Fleet Parameters ###

Fleet<- NULL

Fleet$MinSizeCaught<- LengthAtAge(1.1,Fish,0)

Fleet$MaxSizeCaught<- .95*Fish$Linf