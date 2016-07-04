## Import Libraries into R
library(svDialogs)
library(ggplot2)
library(extrafont)


## Import Data into Data Frame
## Must have column names Pesticide, Dose, Dead, Treated)
File_Name=file.choose()
mainDir=dirname(File_Name)
setwd(mainDir)
Data=read.csv(File_Name,header=T)
Data=Data[complete.cases(Data),]

### Failed (so far) Abbots correction script
{



for (P in Data$Pesticide)
  pest.sub=subset(Data,Pesticide==P)
  #cor.data=matrix(ncol=7,nrow=length(unique(pest.sub$Dose))*length(unique(pest.sub$Genotype)))
  #assign(paste("cor.data",P,sep="."),cor.data)
  for (G in Data$Genotye){
    gen.sub=subset(pest.sub,Genotype==G)
      for (D in Data$Dose){
        dose.sub=subset(gen.sub,Dose==0 | Dose==D)
      sum=ddply(dose.sub,c("Dose"),summarise,
            m=mean(Fraction.Mortality),
            n=length(Fraction.Mortality),
            v=var(Fraction.Mortality))
      t.value=qt(.95,min(sum$n)-1)
      g.value=sum[which(sum$Dose==0),"v"]*(t.value^2)/(((1-sum[which(sum$Dose==0),"m"])^2)*sum[which(sum$Dose==0),"n"])
      p.corr.mort=1-((1-sum[which(sum$Dose==D),"m"])/((1-sum[which(sum$Dose==0),"m"])/(1-g.value)))
      c.i= ((((1-g.value) * (sum[which(sum$Dose==D),"v"]/sum[which(sum$Dose==0),"n"])) + ((((1-sum[which(sum$Dose==D),"m"])^2) * sum[which(sum$Dose==0),"n"])))^.5) * (t.value/((1-sum[which(sum$Dose==0),"m"])/ (1-g)))
      
      
      output[rows,] <- c(j, i, experiment.count, experiment.mean, experiment.var, p.corr.mort, c.i)
   
        
      control.sub=subset(Data,Dose==0 & Genotype==G)
      control.sub.mort=control.sub$Fraction.Mortality
      control.sub.count=length(control.sub.mort)
      control.sub.count=length(control.sub.mort)
      }
  }
}
    
    
## Calculate log dose
Data$modDose=(1000*Data$Dose)+1
conf.level=as.numeric(dlgInput(message="Please Enter What your confidence Interval Should Be (e.g .95 for 95% CI)")$res) 
LD.level=as.numeric(unlist(strsplit(as.character(dlgInput(message="Please Enter What your LD value Should Be (e.g for LD50, write 50")$res),split=",")))
control=as.character(dlgInput(message="Please Enter Which Control Line You Are Using")$res)


GraphDir=paste(getwd(),"/",Sys.Date(),sep="")
dir.create(GraphDir)
setwd(GraphDir)

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

### LD50 function from Paper
LD <- function(r, n, d, conf.level) {
	## Set up a number series 
p=LD.level

## r=number responding, n=number treated, d=dose (untransformed), confidence interval level, 
	mod <- glm(cbind(r, (n-r)) ~ log10(d), family = binomial(link=probit))
	### Calculate heterogeneity correction to confidence intervals according to Finney, 1971, (p.
### 72, eq. 4.27; also called "h")
	het = deviance(mod)/df.residual(mod)
	if(het < 1){het = 1} ### Heterogeneity cannot be less than 1

	## Extract slope and intercept
	summary <- summary(mod, dispersion=het, cor = F)
	intercept <- summary$coefficients[1]
	interceptSE <- summary$coefficients[3]
	slope <- summary$coefficients[2]
	slopeSE <- summary$coefficients[4]
	z.value <- summary$coefficients[6]
	N <- sum(n)
	
	## Intercept and slope (alpha and beta)
	b0<-intercept
	b1<-slope
	
## Slope variance 
	vcov = summary(mod)$cov.unscaled
	var.b0<-vcov[1,1]
	## Intercept variance
	var.b1<-vcov[2,2]
	## Slope intercept covariance
	cov.b0.b1<-vcov[1,2]
	
	## Adjust alpha depending on heterogeneity (Finney, 1971, p. 76)
	alpha<-1-conf.level
	if(het > 1) {talpha <- -qt(alpha/2, df=df.residual(mod))} else {talpha <- -qnorm(alpha/2)}
	
	## Calculate g (Finney, 1971, p 78, eq. 4.36)  
## "With almost all good sets of data, g will be substantially smaller than 1.0 and 
## seldom greater than 0.4."
	g <- het * ((talpha^2 * var.b1)/b1^2)

## Calculate theta.hat for all LD levels based on probits in eta (Robertson et al., 2007, pg. 
## 27; or "m" in Finney, 1971, p. 78)
	eta = family(mod)$linkfun(p/100)  #probit distribution curve
	theta.hat <- (eta - b0)/b1

	## Calculate correction of fiducial limits according to Fieller method (Finney, 1971, 
## p. 78-79. eq. 4.35) 
const1 <- (g/(1-g))*(theta.hat + cov.b0.b1/var.b1) # const1 <- (g/(1-g))*(theta.hat -   cov.b0.b1/var.b1)
const2a <- var.b0 + 2*cov.b0.b1*theta.hat + var.b1*theta.hat^2 - g*(var.b0 - (cov.b0.b1^2/var.b1))
	const2 <- talpha/((1-g)*b1) * sqrt(het * (const2a))
	
	## Calculate the confidence intervals LCL=lower, UCL=upper (Finney, 1971, p. 78-79. eq. 4.35) 
	LCL <- (theta.hat + const1 - const2)
	UCL <- (theta.hat + const1 + const2)
			
	## Calculate variance for theta.hat (Robertson et al., 2007, pg. 27)
	var.theta.hat <- (1/(theta.hat^2)) * ( var.b0 + 2*cov.b0.b1*theta.hat + var.b1*theta.hat^2 )
	
	## Make a data frame from the data at all the different values
	ECtable <- data.frame(
	  "p"=p,
	  "N"=N,
	  "EC"=(10^theta.hat),
	  "LCL"=10^LCL,
	  "UCL"=10^UCL)
	
	## Select output level
	return(ECtable)
}

sd.lc=function(df,LD.level){
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
  
  




## Calculate LC50 values (lc.df) and resistance ratios (RR.df)
{

lc.mat=matrix(nrow=length(unique(Data$Pesticide))*length(unique(Data$Genotype)),ncol=5)
colnames(lc.mat)=c("Pesticide","Genotype","LC50","LCL","UCL")
rownames(lc.mat)=apply(expand.grid(unique(Data$Pesticide), unique(Data$Genotype)), 1, paste, collapse=".")

for (P in unique(Data$Pesticide)){
  for (G in unique(Data$Genotype)){
    Data.sub=subset(Data,Genotype==G & Pesticide==P)
    summary.ld=LD(Data.sub$Dead,Data.sub$Total,Data.sub$modDose,.95)
    lc.mat[paste(P,G,sep="."),]=c(P,G,summary.ld[which(summary.ld[,'p']=="50"),"EC"],summary.ld[which(summary.ld[,'p']=="50"),"LCL"],summary.ld[which(summary.ld[,'p']=="50"),"UCL"])
  }
}
lc.df=data.frame(lc.mat)
lc.df$LC50=as.numeric(as.character(lc.df$LC50))
lc.df$LCL=as.numeric(as.character(lc.df$LCL))
lc.df$UCL=as.numeric(as.character(lc.df$UCL))
write.csv(lc.df,file="Lc50_Summary.csv")
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


}

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






### Individual Dose response in ggplot
for (P in as.character(unique(Data$Pesticide))){

sub.Data=subset(Data,Pesticide==P) 
#pdf(file=paste(p,"Dose_Response.pdf",sep="_"),width=10) 
gp=ggplot(data=sub.Data,aes(x=Dose,y=Survivors,group=interaction(Dose,Genotype)))
gp=gp+geom_dotplot(aes(fill=Genotype),binaxis='y',stackdir="center", binwidth=1,colour="black",dotsize=1,position="dodge")
gp=gp+ggtitle(P) 
gp=gp+ stat_summary(aes(colour=Genotype),fun.data = mean_cl_normal, geom = "crossbar",position=position_dodge()) 
gp=gp+theme(text=element_text(size=18,family="Arial Narrow",face="bold"),
      axis.text.x=element_text(angle = 60, hjust = 1,size=14),
      axis.text.y=element_text(size=14,family="",face="bold"),
      panel.background=element_rect(fill="grey95"),
      panel.border=element_rect(colour="black",fill=NA),
      panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),legend.position="none")
#gp=gp+
print(gp)

#dev.off()
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


