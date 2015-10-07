rm(list=ls())
setwd("C:/Users/lthomas/GitHub/DPSA/Curacao")






sapply(list.files(pattern="[.]R$", path="Functions", full.names=TRUE), source)
source("SubFunctions.R") #Pull in helper functions for assessment modules

Species<-'Acanthurus chirurgus'
Assessment <- 'Doctorfish'
dir.create(Assessment)
Sites<-c("Curacao")
AssessmentName <- paste(Assessment,sep='_')
Directory<- paste(Assessment, "/", sep='')

Font <- 'Helvetica'

FontColor <- 'Black'

PlotFontSize <- 11

Files <- list.files(Assessment)

if ( any(grepl('_LengthData',Files)) )
{
  LengthData <- read.csv(paste(Assessment,'/',Files[grepl('_LengthData',Files)],sep=''), stringsAsFactors = F)
}



source(paste("Doctorfish/Doctorfish_ControlFile.R"))


FigureFolder<- paste(Directory,'Figures/',sep='')

if (file.exists(FigureFolder)==F)
{
  dir.create(FigureFolder,recursive=T)
  
}

Theme<- theme(legend.position='right',plot.background=element_rect(color="white"),
              rect=element_rect(fill='transparent',color=NA),
              text=element_text(size=12,family=Font,color=FontColor),plot.title=element_text(face="italic"),
              axis.text.x=element_text(color=FontColor), axis.text.y=element_text(color="black"),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="black"))

if (exists('LengthData')) {PlotLengthData(LengthData,FigureFolder,Fish,Species,Sites,Theme)}
graphics.off()
