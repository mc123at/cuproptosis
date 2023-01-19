
library("glmnet")
library("survival")
#library("iilasso")
#setwd("D:\\biowolf\\metabolism\\15.lasso")                #设置工作目录
rt=read.table("surSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)            #读取文件

#rt$futime=rt$futime/12
#构建模型1
set.seed(41113256)
x=as.matrix(rt[,c(3:ncol(rt))])
#x=scale(x)
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox")
cvfit=cv.glmnet(x, y, family="cox")

pdf("lambda.pdf", height =5, width = 5) # 图片保存到pdf中
plot(fit, xvar = "lambda", label = T)
dev.off() #关闭图片保存

pdf("min.pdf", height =5, width = 5) #将图片写入到pdf中
plot(cvfit)
dev.off() # 关闭图片，图片保存到pdf中


#plot(fit)
#plot(fit, xvar = "lambda", label = TRUE)
#plot(cvfit)

#输出相关基因系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)


#输出train组风险值
#rt$futime=rt$futime/365
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.1.LUAD_Risk.txt",sep="\t",quote=F,row.names=F)


#输出tain 5 yesr  组风险值
rt=read.table("2.valid.geoExpTime.txt", sep = "\t",header = T, na.strings = "NA",stringsAsFactors=FALSE, check.names=F,row.names=1,quote="") #基因名有'引号的需要加上quote=""
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
#rt$futime=rt$futime/365
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
#risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.2.GSE72094.txt",sep="\t",quote=F,row.names=F)



#输出GSE72094组风险值
rt=read.table("3.valid.geoExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
#rt$futime=rt$futime/365
test2=rt[,lassoGene]
testScore2=apply(test2,1,myFun)
outCol=c("futime","fustat",lassoGene)

risk=as.vector(ifelse(testScore2>median(testScore2),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore2),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.3.GSE68465.txt",sep="\t",quote=F,row.names=F)



















#输出GSE72094 5yesrs组风险值
rt=read.table("4.expPFI.5years.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
rt$futime=rt$futime/365
train5=rt[,lassoGene]
train5score=apply(train5,1,myFun)
outCol=c("futime","fustat",lassoGene)

risk=as.vector(ifelse(train5score>median(train5score),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(train5score),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.4.train.expPFI.5years.txt",sep="\t",quote=F,row.names=F)



#输出GSE68465 组风险值
rt=read.table("5.GSE72094.ExpOs.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
rt$futime=rt$futime/365
test5=rt[,lassoGene]
test5score=apply(test5,1,myFun)
outCol=c("futime","fustat",lassoGene)

risk=as.vector(ifelse(test5score>median(test5score),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(test5score),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.5.GSE72094.ExpOs.txt",sep="\t",quote=F,row.names=F)




#输出GSE68465 5years 组风险值
rt=read.table("6.input.GSE72094.5years.ExpOs.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
rt$futime=rt$futime/365
vio1=rt[,lassoGene]
vio1score=apply(vio1,1,myFun)
outCol=c("futime","fustat",lassoGene)

risk=as.vector(ifelse(vio1score>median(vio1score),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(vio1score),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.6.GSE72094.ExpOs.5years.txt",sep="\t",quote=F,row.names=F)





#输出tcgapfi 组风险值
rt=read.table("7.GSE68465ExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
#rt$futime=rt$futime/1
GSE30589gene=rt[,lassoGene]
GSE30589score=apply(GSE30589gene,1,myFun)
outCol=c("futime","fustat",lassoGene)

risk=as.vector(ifelse(GSE30589score>median(GSE30589score),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(GSE30589score),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.7.GSE68465ExpTime.txt",sep="\t",quote=F,row.names=F)





#输出tcgapfi 5 yrs组风险值
rt=read.table("8.GSE68465.OS.5years.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
#rt$futime=rt$futime/365
wholeFinalGeneExp=rt[,lassoGene]
wholeScore=apply(wholeFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)

risk=as.vector(ifelse(wholeScore>median(wholeScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(wholeScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="RISK.8.GSE68465.OS.5years.txt",sep="\t",quote=F,row.names=F)










#输出immune535组风险值
rt=read.table("ExpTime.Immune.36.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
rt$futime=rt$futime/365
immune535GeneExp=rt[,lassoGene]
immune535Score=apply(immune535GeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)

risk=as.vector(ifelse(immune535Score>median(immune535Score),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(immune535Score),risk)
write.table(cbind(id=rownames(outTab),outTab),file="immune36Risk.txt",sep="\t",quote=F,row.names=F)



