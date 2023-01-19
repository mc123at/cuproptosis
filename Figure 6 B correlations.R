#环境设置
library(reshape2)
library(ggplot2)
library(cowplot)
## 
## ********************************************************
## Note: As of version 1.0.0, cowplot does not change the
##   default ggplot2 theme anymore. To recover the previous
##   behavior, execute:
##   theme_set(theme_cowplot())
## ********************************************************
library(plyr)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#参数设置
symbol <- "Riskscore"
topnumber = 22  #取正相关排名Top X和负相关排名Top X的基因

# plot with special pointnumber by groupping gene expression levels into equally sized bins
pointnumber = 50 
# plot without bins
#pointnumber = NULL

#输入文件
#tcga_expr <- read.csv("TMB and GNPNAT1.txt", row.names = 1)

tcga_expr=read.table("CIBERSORT-Results.downFromGDC.txt", header=T,sep="\t",row.names=1,check.names=F)
tcga_expr[1:3,1:3]

#计算相关系数
# 目标基因
target.exps <- tcga_expr[symbol,]
# 其余基因
other.expr <- tcga_expr[-which(rownames(tcga_expr)==symbol),]



# pearson
sgcor <- as.data.frame(cor(t(other.expr), t(target.exps))) #默认pearson
colnames(sgcor) <- "r_pearson"
# 计算pvalue，运行时间较长
sgcor$pval_pearson <- apply(other.expr, 1, function(x) (cor.test(x, t(target.exps))$p.value))



# spearman
sgcor_spearman <- as.data.frame(cor(t(other.expr), t(target.exps), method = "spearman"))
colnames(sgcor_spearman) <- "r_spearman"
sgcor_spearman$pval_spearman <- apply(other.expr, 1, function(x)(cor.test(x, t(target.exps), method = "spearman")$p.value))






# 把相关系数、pvalue都写进cors里
cors <- cbind(sgcor, sgcor_spearman)
cors$gene <- rownames(other.expr)
head(cors)
write.table(cors,file="1.allCOR.xls",sep="\t",col.names = NA, row.names=T,quote=F) #第一列列名为空，这里NA可以解决




# 取相关系数排名前几位的基因，在参数设置里修改topnumber
# 此处取前X和后X，也就是正相关TopX和负相关TopX
newcor <- cors[!(is.na(cors$r_pearson)),]
dim(newcor)

sortcor <- newcor[order(newcor$r_pearson, newcor$r_spearman, decreasing = T),]
topcor <- sortcor[c(1:topnumber, #正相关Top X
                    (nrow(sortcor) - topnumber + 1):nrow(sortcor)),] #负相关Top X
rownames(topcor)

# 提取相关系数排名前几位基因的表达矩阵,目标基因也会写入table中为第一列基因
genes <- c(symbol,rownames(topcor))
genesexps <- as.data.frame(t(tcga_expr[genes,]))
sortgenesexps <- genesexps[order(genesexps[,1]),]
write.table(sortgenesexps,file="2.topGeneExp.xls",sep="\t",col.names = NA, row.names=T,quote=F) #第一列列名为空，这里NA可以解决


samplenum <- nrow(sortgenesexps)



if(is.null(pointnumber)){
  pointnumber=samplenum
}

# plot with special pointnumber by groupping gene expression levels into equally sized bins
group <- as.integer(cut(1:samplenum, breaks=c(seq(from=0.5, to=samplenum, by=samplenum/pointnumber), samplenum+0.5)))
ddf <- data.frame(row.names = 1:pointnumber)

for( i in 1:(1 + topnumber*2)){
  ddf <- cbind(ddf,tapply(sortgenesexps[,i],group,median))
}

colnames(ddf) <- c(symbol,topcor$gene)
mddf <- melt(ddf,id.vars=symbol)

# 在图中显示pearson的r和P
# 或者显示spearman的r和P，就换成r_spearman和pval_spearman
mddf$r <- topcor[mddf$variable,]$r_pearson
mddf$P <- topcor[mddf$variable,]$pval_pearson 

#mddf$r <- topcor[mddf$variable,]$r_spearman
#mddf$P <- topcor[mddf$variable,]$pval_spearman 










#开始画图

#批量画图

plist <- dlply(mddf, .(variable), function(trig){ggplot(trig, aes_string(x=symbol, y="value")) +
    geom_point() +
    ylab(unique(trig$variable)) +
    theme_classic() +
    ggtitle(paste0("r = ", round(unique(trig$r),2),
                   "\nP = ", sprintf("%1.1e", unique(trig$P)))) +
    geom_smooth(method = "lm", se=T, colour = "#0086b3")})

pg <- plot_grid(plotlist = plist, ncol=5, align = "hv")
pg
filename <- paste0(symbol, "_cor.pdf")
ggsave(filename, width =9, height = 9)






