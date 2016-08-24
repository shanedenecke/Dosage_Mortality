



##############################################################################################################################################


dir.create(file.path(getwd(),"Resistace_Ratio_Graphs"))
setwd(file.path(getwd(),"Resistace_Ratio_Graphs"))
RR.df=data.frame()
lc.df=arrange(lc.df,Pesticide)
for (P in unique(lc.df$Pesticide)){
    lc.sub=subset(lc.df, Pesticide==P)  
    for (G in unique(lc.df$Genotype)){  
        RR=lc.sub$LC50[which(lc.sub$Genotype==G)]/lc.sub$LC50[which(lc.sub$Genotype==control)]
        RR.LCL=RR*min(lc.sub$LCL[which(lc.sub$Genotype==G)]/lc.sub$LC50[which(lc.sub$Genotype==G)],lc.sub$LCL[which(lc.sub$Genotype==control)]/lc.sub$LC50[which(lc.sub$Genotype==control)])
        RR.UCL=RR*max(lc.sub$UCL[which(lc.sub$Genotype==G)]/lc.sub$LC50[which(lc.sub$Genotype==G)],lc.sub$UCL[which(lc.sub$Genotype==control)]/lc.sub$LC50[which(lc.sub$Genotype==control)])
        RRs=c(RR,RR.LCL,RR.UCL) 
        RR.df=rbind(RR.df,RRs)
    }} 
colnames(RR.df)=c("RR","RR.LCL","RR.UCL")

sum.df=cbind(lc.df,RR.df)





### Plot Resistance Ratios for All insecticides. Not quite working yet
{
    
    RR.plot.df=arrange(sum.df,RR)
    RR.plot.df$Pesticide=factor(RR.plot.df$Pesticide,levels=unique(RR.plot.df$Pesticide))
    RR.plot.color=colorRampPalette(c("orangered1","steelblue1"))(length(RR.plot.df$Pesticide))
    
    pdf(file="Resistance_Ratio_Graph.pdf")
    RRplot=ggplot(RR.plot.df,aes(x=Pesticide,y=RR,fill=Genotype))
    RRplot=RRplot+geom_bar(stat="identity",position="dodge")
    RRplot=RRplot+geom_errorbar(aes(ymin=RR.LCL,ymax=RR.UCL),position="dodge",stat="identity")
    RRplot=RRplot+ylim(c(0,RR.UCL*1.5))
    RRplot=RRplot+theme(text=element_text(size=18,face="bold"),
                        axis.text.x=element_text(angle = 60, hjust = 1,size=14),
                        axis.text.y=element_text(size=14,family="",face="bold"),
                        panel.background=element_rect(fill="grey95"),
                        panel.border=element_rect(colour="black",fill=NA),
                        panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())
    RRplot=RRplot+geom_hline(yintercept=1,linetype=2)
    print(RRplot)
    dev.off()
}


##LC50 Comparisons in ggplot2
{
    
    for (P in unique(lc.df$Pesticide)){ 
        
        sub.lc=subset(lc.df,Pesticide==P)
        
        #pdf(file=paste(p,"LC_50_Comparison.pdf",sep="_"),width=10) 
        lc.plot=ggplot(data=sub.lc,aes(x=Genotype,y=LC50,fill=Genotype)) 
        lc.plot=lc.plot+geom_bar(stat="identity",position="dodge")
        lc.plot=lc.plot+geom_errorbar(aes(ymin=LCL,ymax=UCL))
        lc.plot=lc.plot+ylim(c(0,max(as.numeric(sub.lc$LC50)*1.5)))
        lc.plot=lc.plot+theme(text=element_text(size=18,family="Arial Narrow",face="bold"),
                              axis.text.x=element_text(angle = 60, hjust = 1,size=14),
                              axis.text.y=element_text(size=14,family="",face="bold"),
                              panel.background=element_rect(fill="grey95"),
                              panel.border=element_rect(colour="black",fill=NA),
                              panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),legend.position="none")
        lc.plot=lc.plot+scale_fill_manual(values=c("orange","black"))
        print(lc.plot)
    }
}





### Base Plot Graphs for LC50 comparisons
{{
    attach(sub.lc,warn.conflicts=FALSE)
    pdf(file=paste(P,"LC50_barplot.pdf",sep="_"),width=10)
    
    par(mar=c(2,2,2,2),mgp=c(3.2,1,0),mfrow=c(1,1))
    bp=barplot(LC50,col=c("blue","red"),ylim=c(0,max(LC50)+.4*max(LC50)),family="serif", 
               ylab=list("LC50 Value (ppb)",font=2,cex=1.8),main=paste(P),cex.main=2,axes=F,las=2)
    segments(bp,LC50, bp,UCL,col="black",lwd=2)
    segments(x0=bp-.3, y0=UCL,x1=bp+.3,lwd=2,col="black")
    segments(bp,LC50, bp,UCL,col="black",lwd=2)
    segments(x0=bp-.3, y0=LCL,x1=bp+.3,lwd=2,col="black")
    axis(side=2,at=seq(0,signif(max(UCL),1),length.out=6),font=2,cex=1.2,las=1,family="serif")
    text(x=bp-.1,y=par("usr")[3]-2,srt = 0, xpd = TRUE,adj=0,labels =paste(unique(Genotype)), font=2, family="serif",cex=1.5)
    mtext(side=1,at=median(bp)-.1, text="Genotype",font=2,cex=2.5,line=3,family="serif")
    
    dev.off()
    detach(sub.lc)
}}


}






##sd.lc=function(df,LD.level){
prob=glm(cbind(Dead,Total-Dead)~log10(modDose),family=binomial(link="probit"),data=df)
prob.sum=summary(prob)
dist=family(prob)$linkfun(LD.level/100)
b0=prob.sum$coefficients[1]
b1=prob.sum$coefficients[2]
theta.hat=10^(dist-b0)/b1
lc.table=((10^(dist-b0)/b1)-1)/1000
ll.table=10^dist-confint(prob)[1,1]/confint(prob)[2,1]-1/1000


return(lc.table)
} 


prob=glm(cbind(Dead,Total-Dead)~log10(modDose)+Genotype,family=binomial(link="probit"),data=sub.Pesticide)
wald.test(b=coef(prob), Sigma=vcov(prob), Terms=4:5)

sub.Data=subset(Data,Genotype==G & Pesticide==P)
prob=glm(cbind(Dead,Total-Dead)~log10(modDose),family=binomial(link="probit"),data=sub.Data)
prob.sum=summary(prob,dispersion=het,cor=F)
p=seq(1,99,1)
eta = family(prob)$linkfun(p/100)  #probit distribution curve
b0=prob.sum$coefficients[1]
b1=prob.sum$coefficients[2]
theta.hat=(eta - b0)/b1
l50=((10^(dist-b0)/b1)/1000)-1
}
