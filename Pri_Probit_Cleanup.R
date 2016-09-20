
File_Name <- file.choose()
mainDir <- dirname(File_Name)
setwd(mainDir)
raw.Data <- read.csv(File_Name,header=T,stringsAsFactors = F)
raw.Data <- raw.Data[complete.cases(raw.Data),]
GraphDir=paste(getwd(),"/",format(Sys.Date(),format="%d %b %y"),"_Pri",sep="")
dir.create(GraphDir)
setwd(GraphDir)



for (G in unique(raw.Data$Genotype)){
    for (P in unique(raw.Data$Pesticide)){
        sub.data <- subset(raw.Data,Pesticide==P & Genotype==G)
        
        pri <- cbind(sub.data$Dose,sub.data$Total,sub.data$Dead)
        
        write.table(pri,file=paste("pri",P,G,".txt",sep="_"),sep="\t",col.names=F,row.names=F)
    }}