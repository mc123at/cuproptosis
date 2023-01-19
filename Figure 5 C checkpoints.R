
#引用包
library(reshape2)
library(ggpubr)
#expFile="LUAD.Tumor.Exp.txt"                 #表达输入文件
#geneFile="gene.txt"                 #基因列表文件
#scoreFile="ScoreGroup.txt"      #score文件
#setwd("D:\\biowolf\\ICI\\30.scoreGene")     #设置工作目录

#读取输入文件
#exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
#score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
#gene=read.table(geneFile, header=F, sep="\t", check.names=F)
rt=read.table("sel.immune.LUAD.txt", header=T, sep="\t", check.names=F, row.names=1)
#合并数据
#sameGene=intersect(row.names(exp), as.vector(gene[,1]))
#geneExp=as.data.frame(t(exp[sameGene,]))
#rownames(geneExp)=gsub("(.*?)\\_(.*?)", "\\2", rownames(geneExp))
#sameSample=intersect(rownames(geneExp), row.names(score))
#geneExp=geneExp[sameSample,]
#score=score[sameSample,]
#rt=cbind(geneExp, ICIscore=score[,"group"])

#把数据转换成ggplot2输入文件
data=melt(rt, id.vars=c("ICIscore"))
colnames(data)=c("ICIscore", "Gene", "Expression")
data$ICIscore=factor(data$ICIscore, levels=c("low", "high"))

#绘制箱线图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(rt[,"ICIscore"])))]
p=ggboxplot(data, x="Gene", y="Expression", color="ICIscore", 
            ylab="Gene expression",
            xlab="",
            legend.title="Risk score",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="boxplot.pdf",width=9,height=4)                          #输出图片文件
p+stat_compare_means(aes(group=ICIscore),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

