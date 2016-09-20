## Import Libraries into R
library(svDialogs)
library(ggplot2)
library(extrafont)
loadfonts()

## Import Data into Data Frame
## Must be formated with columns Pesticide, Genotype rep, Dose, Alive, Dead, Total, 
## Optional column names Life_Stage, Date, Conditions
## Example formatted file can be found by unhashtaging the next two lines
##write.csv(url("******"),file="example_formatted_file.csv")
##head(read.csv("example_formatted_file.csv")

File_Name <- file.choose()
mainDir <- dirname(File_Name)
setwd(mainDir)
raw.Data <- read.csv(File_Name,header=T)
raw.Data <- raw.Data[complete.cases(raw.Data),]
GraphDir=paste(getwd(),"/",format(Sys.Date(),format="%d %b %y"),sep="")
dir.create(GraphDir)
setwd(GraphDir)
 

### Add Control treatments to each insecticide
Data <- data.frame()
for (P in as.character(unique(raw.Data$Pesticide))){
    
    sub.Data.Pesticide=subset(raw.Data,Pesticide==P) 
    
    if(P=="Spinosad" | P=="Imidacloprid" | P=="Chlorantraniliprole" | P=="Nitenpyram" | P=="Malathion"){
        sub.Data.Pesticide <- rbind(sub.Data.Pesticide,subset(raw.Data,Pesticide=="H20"))
        
    } else if (P=="Ivermectin" | P=="Lufenuron" | P=="Sulfoxaflor" | P=="Nit_ver(100)"){
        sub.Data.Pesticide <- rbind(sub.Data.Pesticide,subset(raw.Data,Pesticide=="DMSO"))
        
    } else if (P=="Pyripole"){
        sub.Data.Pesticide <- rbind(sub.Data.Pesticide,subset(raw.Data,Pesticide=="DMSO"))
    } else if (P=="H20" | P=="DMSO"){
        break
    }
    
    sub.Data.Pesticide$Pesticide <- rep(sub.Data.Pesticide$Pesticide[1],length(sub.Data.Pesticide$Pesticide))
    Data <- rbind(sub.Data.Pesticide,Data)
}
    



## LD50 prep stuff and function from Paper ###

    ## Prep data frame and enter variables
{
Data$modDose=(1000*Data$Dose)+1
conf.level=.95#as.numeric(dlgInput(message="Please Enter What your confidence Interval Should Be (e.g .95 for 95% CI)")$res) 
LD.level=50#as.numeric(unlist(strsplit(as.character(dlgInput(message="Please Enter What your LD value Should Be (e.g for LD50, write 50")$res),split=",")))
nm=as.character(dlgInput(message="Please Enter A unique name")$res)
}

    ## LD50 function from paper
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

    ## Calculate LC50 values for each genotype and pesticide
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
}  
lc.df[,3:5] <- (lc.df[,3:5]-1)/1000



## Individual Dose Dotplot
dir.create(file.path(getwd(),"Indiviual_Vial_Dotplots"))
setwd(file.path(getwd(),"Indiviual_Vial_Dotplots"))
for (P in as.character(unique(Data$Pesticide))){
    
    sub.Data=subset(Data,Pesticide==P) 
    
    sub.Data$Dose <- as.factor(as.character(sub.Data$Dose))
    
    pdf(file=paste(P,nm,"Dose_Response.pdf",sep="_"),width=10) 
    gp=ggplot(data=sub.Data,aes(x=Dose,y=Alive,group=interaction(Dose,Genotype)))
    gp=gp+geom_dotplot(aes(fill=Genotype),binaxis='y',stackdir="center", binwidth=.1,colour="black",position=position_dodge(width=.5),dotsize = 10)
    gp=gp+ggtitle(paste("Individual Vial Response for",P,sep=" ")) 
    gp=gp+stat_summary(geom="errorbar",colour='black',size=.6,fun.data=mean_cl_normal,position=position_dodge(width=.5),width=.5)
    gp=gp+ylab("Number Emerging\n")
    gp=gp+xlab("\nDose (ppm)")
    gp=gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar",colour='black',position=position_dodge(width=.5),width=.5)
    gp=gp+scale_y_continuous(limits=c(0,max(as.numeric(sub.Data$Alive))))
    gp=gp+scale_x_discrete(breaks=as.factor(sub.Data$Dose))
    
    gp=gp+theme(text=element_text(size=18 ,face="bold"),
                axis.text.x=element_text(angle = 60, hjust = 1,size=14),
                axis.text.y=element_text(size=14,family="",face="bold"),
                panel.background=element_rect(fill="grey95"),
                panel.border=element_rect(colour="black",fill=NA),
                panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())
    print(gp)
    
    dev.off()
} 
setwd(GraphDir)

### Abbots Correction
dir.create(file.path(getwd(),"Abbots Correction"))
setwd(file.path(getwd(),"Abbots Correction"))
for (P in as.character(unique(Data$Pesticide))){
    
    sub.Data.Pesticide=subset(Data,Pesticide==P) 
    w <- max(sub.Data.Pesticide$Dose)/5
    sum.data <- matrix(ncol=4,nrow=1)
    colnames(sum.data) <- c("Genotype","Dose","Corrected.Mortality","CI")
    
    for (G in unique(sub.Data.Pesticide$Genotype)){
         sub.Data.Genotype <- subset(sub.Data.Pesticide,Genotype==G)
        for (D in unique(sub.Data.Pesticide$Dose)){
            
            ## Calculate Abbots Correction
            control.mort <- sub.Data.Genotype$Percent_Mortality[sub.Data.Genotype$Dose==0]
            control.count <- length(control.mort)
            control.mean <- mean(control.mort)
            control.var <- var(control.mort)
            ex.mort <- sub.Data.Genotype$Percent_Mortality[sub.Data.Genotype$Dose==D]
            ex.count <- length(ex.mort)
            ex.mean <- mean(ex.mort)
            ex.var <- var(ex.mort)
            min.samples <- min(control.count,ex.count)
            dof <- min.samples - 1 #degrees of freedom
            t.value <- qt(.05, dof)
            g <- (control.var * (t.value^2))/(((1-control.mean)^2) * control.count)
            p.corr.mort <- 1 - ((1-ex.mean)/((1-control.mean)/(1-g)))
            c.i <- ((((1-g) * (ex.var/ex.count)) + ((((1-ex.mean)^2) * control.var) / (((1-control.mean)^2) * control.count)))^0.5) * (t.value/((1-control.mean) / (1-g)))
            temp.vec <- c(G,D,p.corr.mort,c.i)
            sum.data <- rbind(temp.vec,sum.data)
        }}   
          rownames(sum.data) <- NULL
          sum.data <- data.frame(sum.data[complete.cases(sum.data),],stringsAsFactors = FALSE)
          sum.data[,2:4] <- as.numeric(as.character(unlist(sum.data[,2:4])))
          sum.data[["Corrected.Survival"]] <- (1-sum.data$Corrected.Mortality)*100
          sum.data[["CI"]] <-sum.data[["CI"]]*100 
          write.csv(sum.data,file=paste(P,"Corrected_Data.csv",sep="_"))
        sum.data$Dose <- as.factor(as.character(sum.data$Dose)) 
            ## Plot Abbots Correction 
        pdf(file=paste(P,nm,"Abbots_Correction.pdf",sep="_"),width=10) 
        gp=ggplot(data=sum.data,aes(x=Dose,y=Corrected.Survival,group=interaction(Dose,Genotype)))
        gp=gp+geom_errorbar(aes(ymin=Corrected.Survival,ymax=Corrected.Survival-CI),position=position_dodge(width=.5),stat="identity",width=.5)   
        gp=gp+geom_bar(aes(fill=Genotype),position=position_dodge(width=.5),stat="identity",colour="black",width=.5)
        gp=gp+ggtitle(paste("Abbots Correction",P,sep=" "))  
        gp=gp+scale_y_continuous(breaks=seq(0,100,by=20),limits=c(0,140))
        gp=gp+ylab("Corrected Percent Survival\n")
        gp=gp+scale_x_discrete(breaks=sum.data$Dose)
        gp=gp+xlab("\nDose (ppm)")
        gp=gp+scale_x_discrete(breaks=unique(sum.data$Dose))
        gp=gp+theme(text=element_text(size=18,face="bold"),
        axis.text.x=element_text(angle = 60, hjust = 1,size=14),
        axis.text.y=element_text(size=14,family="",face="bold"),
        panel.background=element_rect(fill="grey95"),
        panel.border=element_rect(colour="black",fill=NA),
        panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())
                    print(gp) 
                    
                    dev.off()
}  
setwd(GraphDir)        


## Plot Lc50 Graphs for each insecticide
dir.create(file.path(getwd(),"Lc50_Comparison"))
setwd(file.path(getwd(),"Lc50_Comparison"))
write.csv(lc.df,file="Lc50_Summary.csv")
for (P in unique(lc.df$Pesticide)){ 
    sub.lc=subset(lc.df,Pesticide==P)
    limits <- c(0,max(as.numeric(sub.lc$UCL)*1.1))
    pdf(file=paste(P,"LC_50_Comparison.pdf",sep="_"),width=10) 
    lc.plot=ggplot(data=sub.lc,aes(x=Genotype,y=LC50,fill=Genotype)) 
    lc.plot=lc.plot+geom_bar(stat="identity",position="dodge",width=.4)
    lc.plot=lc.plot+geom_errorbar(aes(ymin=LCL,ymax=UCL),width=.25)
    lc.plot=lc.plot+ylim(limits)
    lc.plot=lc.plot+ggtitle(paste("LC50 Values for",P,sep=" "))
    lc.plot=lc.plot+ylab("Lc50 Value (ppm)\n")
    lc.plot=lc.plot+xlab("\nGenotype")
    lc.plot=lc.plot+theme(text=element_text(size=18,face="bold"),
                          axis.text.x=element_text(angle = 60, hjust = 1,size=14),
                          axis.text.y=element_text(size=14,family=""),
                          panel.background=element_rect(fill="grey95"),
                          panel.border=element_rect(colour="black",fill=NA),
                          panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),legend.position="none")
    print(lc.plot)
dev.off()
}
setwd(GraphDir)



