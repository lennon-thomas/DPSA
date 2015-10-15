library(ggplot2)
  
setwd("C:/Users/lthomas/GitHub/DPSA/Montserrat")

data<-read.csv("results.csv")

fplot<-ggplot(data=data,aes(x=Species, y=F,fill=Method))+
              geom_bar(stat="identity",position=position_dodge(),color="black")+
  geom_hline(yintercept=mean(1),linetype="dashed",size=1.0)+
labs(x="",y=expression(F/F[MSY]))+
  #geom_errorbar(aes(ymin=data$lower, ymax=data$upper), width=.2,position=position_dodge(.9),colour="black")+
  scale_y_continuous(breaks=seq(0,12,1), expand = c(0, 0))
                                                                                  
fplot<-fplot+theme(
  axis.text = element_text(size = 14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white"),
  axis.line=element_line(colour="Black"),axis.text.x = element_text(angle = 45,size=14,hjust=1),axis.title.y=element_text(size=16,vjust=1.5,face="italic"),
  axis.title.x=element_text(size=16,vjust=0.1),legend.position=c(0.87,0.9),legend.title=element_text(size=14),legend.text=element_text(size=12))

fplot<-fplot + 
  
  scale_fill_manual(values=c("#666666", "#CCCCCC", "white"), #"#0066CC", "#FF9999", "#00CC66"
                       name="Assessment Method",
                       breaks=c("Catch MSY", "Lbar", "lbspr"),
                       labels=c("Catch MSY", "LBAR", "LB-SPR"))


plot(fplot)

pdf(file="Fmsy_bw.pdf")
fplot
dev.off()
