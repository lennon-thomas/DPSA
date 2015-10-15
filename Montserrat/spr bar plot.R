

rm(list = ls())
library(ggplot2)


setwd("C:/Users/lthomas/GitHub/DPSA/Montserrat")
data<-read.csv("SPR results.csv")
names(data)
data$Standard<-rep.int(30,6)
Font <- 'Helvetica'

FontColor <- 'Black'

PlotFontSize <- 11

results_spr<-ggplot(data=data, aes(x=Species.Name, y=SPR,fill=as.factor(Type)))+
 
  stat_summary(fun.y=mean, geom="bar",colour="black",aes(width=0.5))+scale_fill_manual(values=c("#666666", "#CCCCCC", "white"))+
  geom_errorbar(aes(ymin=data$Lower, ymax=data$Upper), width=.2,
                position=position_dodge(.9)) +                
  geom_hline(yintercept=mean(data$Standard),linetype="dashed")+
    labs(x="", y="Spawning Potential Ratio (SPR) %")+                 #scale_x_continuous(limits=c(1986,2014), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,100), expand = c(0, 0))

results_spr<-results_spr+        theme(
          axis.text = element_text(size = 14),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line=element_line(colour="Black"),axis.text.x = element_text(angle = 45,size=14,hjust=1),axis.title.y=element_text(size=16,vjust=1.5),
          axis.title.x=element_text(size=16,vjust=0.1),legend.position="none"
)

plot(results_spr)

pdf(file="SPR Results.pdf")
results_spr
dev.off()
