library(reshape2)
library(ggpubr)
rt=read.table("sel.immune.LUAD.txt", header=T, sep="\t", check.names=F, row.names=1)
data=melt(rt, id.vars=c("ICIscore"))
colnames(data)=c("ICIscore", "Gene", "Expression")
data$ICIscore=factor(data$ICIscore, levels=c("low", "high"))
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(rt[,"ICIscore"])))]
p=ggboxplot(data, x="Gene", y="Expression", color="ICIscore", 
            ylab="Gene expression",
            xlab="",
            legend.title="Risk score",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="boxplot.pdf",width=9,height=4)              
p+stat_compare_means(aes(group=ICIscore),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()