## Import Libraries into R
library(svDialogs)
library(ggplot2)


## Import Data into Data Frame
## Must have column names Pesticide, Dose, Dead, Treated)
File_Name=file.choose()
mainDir=dirname(File_Name)
setwd(mainDir)
Data=read.csv(File_Name,header=T)
Data=Data[complete.cases(Data),]


## Calculate log dose
Data$modDose=(1000*Data$Dose)+1
conf.level=as.numeric(dlgInput(message="Please Enter What your confidence Interval Should Be (e.g .95 for 95% CI)")$res) 
LD.level=as.numeric(dlgInput(message="Please Enter What your LD value Should Be (e.g for LD50, write 50")$res)

#a=Sys.Date
GraphDir=paste(getwd(),"/",Sys.Date(),sep="")
dir.create(GraphDir)
setwd(GraphDir)

sub=subset(Data,Pesticide=="Imidacloprid")
par(mar=c(2,2,2,2),mgp=c(3,1,0),mfrow=c(1,1),xpd=TRUE)
plot(sub$Survivors~sub$modDose)

### LD50 function
LD <- function(r, n, d, conf.level) {
	## Set up a number series 
#p <- seq(1, 99, 1)

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
	
	## Intercept (alpha)
	b0<-intercept
	## Slope (beta)
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
	  "EC"=10^theta.hat,
	  "LCL"=10^LCL,
	  "UCL"=10^UCL, 
	  "slope"=slope, 
	  "slopeSE"=slopeSE, 
	  "intercept"=intercept, 
	  "interceptSE"=interceptSE, 
	  "z.value"=z.value, 
	  "chisquare"=deviance(mod), 
	  "df"=df.residual(mod), 
	  "h"=het, 
	  "g"=g,
	  "theta.hat"=theta.hat,
	  "var.theta.hat"=var.theta.hat)
	
	## Select output level
	return(ECtable)
}




{
emat=matrix(nrow=length(unique(Data$Pesticide))*length(unique(Data$Genotype)),ncol=5)
colnames(emat)=c("Pesticide","Genotype","LC50","LCL","UCL")
rownames(emat)=apply(expand.grid(unique(Data$Pesticide), unique(Data$Genotype)), 1, paste, collapse=".")

for (p in unique(Data$Pesticide)){
  for (g in unique(Data$Genotype)){
    Data.sub=subset(Data,Genotype==g & Pesticide==p)
    summary=LD(Data.sub$Dead,Data.sub$Total,Data.sub$modDose,.95)
    emat[paste(p,g,sep="."),]=c(p,g,summary[50,"EC"],summary[50,"LCL"],summary[50,"UCL"])
    write.csv(summary,file=paste(p,g,"Summary_Table.csv",sep="_"))
  }
}

pt=data.frame(emat)
pt$LC50=as.numeric(as.character(pt$LC50))
pt$LCL=as.numeric(as.character(pt$LCL))
pt$UCL=as.numeric(as.character(pt$UCL))

write.csv(pt,file="Overall_LC50_Summary.csv")
for (p in unique(emat[,"Pesticide"])){ 
sub.pt=subset(pt,Pesticide==p)
attach(sub.pt)
pdf(file=paste(p,"LC50_barplot.pdf",sep="_"),width=10)

par(mar=c(6,6,3,3),mgp=c(3.2,1,0),mfrow=c(1,1))
bp=barplot(LC50,col=c("blue","red"),ylim=c(0,3+max(UCL)),family="serif", 
ylab=list("LC50 Value (ppb)",font=2,cex=1.8),main=paste(p),cex.main=2,axes=F,las=2)
segments(bp,LC50, bp,UCL,col="black",lwd=2)
segments(x0=bp-.3, y0=UCL,x1=bp+.3,lwd=2,col="black")
segments(bp,LC50, bp,UCL,col="black",lwd=2)
segments(x0=bp-.3, y0=LCL,x1=bp+.3,lwd=2,col="black")
axis(side=2,at=seq(0,signif(max(UCL),1),length.out=6),font=2,cex=1.2,las=1,family="serif")
text(x=bp-.1,y=par("usr")[3]-2,srt = 0, xpd = TRUE,adj=0,labels =paste(unique(Genotype)), font=2, family="serif",cex=1.5)
mtext(side=1,at=median(bp)-.1, text="Genotype",font=2,cex=2.5,line=3,family="serif")

dev.off()
}


for (p in as.character(unique(Data$Pesticide))){

sub.Data=subset(Data,Pesticide==p) 
pdf(file=paste(p,"Dose_Response.pdf",sep="_"),width=10) 
gp=ggplot(data=sub.Data,aes(Dose,Survivors,fill=Genotype)) 
#gp=gp+geom_bar(stat="identity", position = "dodge")
gp=gp+ggtitle(p) 
gp=gp+ stat_summary(fun.data = mean_cl_normal, geom = "errorbar",position=position_dodge()) 
gp=gp+stat_summary(fun.y=mean,position=position_dodge(),geom="bar") 
print(gp)

dev.off()
} 

}



p <- ggplot(mtcars, aes(wt, mpg))
p <- p + geom_point()
print(p)



