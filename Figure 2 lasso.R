library("glmnet")
library("survival")
#library("iilasso")
#setwd("D:\\biowolf\\metabolism\\15.lasso")               
rt=read.table("surSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)           
#rt$futime=rt$futime/12
set.seed(41113256)
x=as.matrix(rt[,c(3:ncol(rt))])
#x=scale(x)
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox")
cvfit=cv.glmnet(x, y, family="cox")
pdf("lambda.pdf", height =5, width = 5) 
plot(fit, xvar = "lambda", label = T)
dev.off() 
pdf("min.pdf", height =5, width = 5) 
plot(cvfit)
dev.off() 
#plot(fit)
#plot(fit, xvar = "lambda", label = TRUE)
#plot(cvfit)
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)
#rt$futime=rt$futime/365
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.1.LUAD_Risk.txt",sep="\t",quote=F,row.names=F)
rt=read.table("2.valid.geoExpTime.txt", sep = "\t",header = T, na.strings = "NA",stringsAsFactors=FALSE, check.names=F,row.names=1,quote="")
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
#rt$futime=rt$futime/365
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
#risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.2.GSE72094.txt",sep="\t",quote=F,row.names=F)
rt=read.table("3.valid.geoExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
#rt$futime=rt$futime/365
test2=rt[,lassoGene]
testScore2=apply(test2,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore2>median(testScore2),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore2),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.3.GSE68465.txt",sep="\t",quote=F,row.names=F)
