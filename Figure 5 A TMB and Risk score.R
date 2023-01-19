library(reshape2)
library(ggplot2)
library(cowplot)
library(plyr)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
symbol <- "riskscore"
topnumber = 1 
pointnumber = NULL
tcga_expr=read.table("RISK.1.Train_OS_Risk.txt", header=T,sep="\t",row.names=1,check.names=F)
tcga_expr[1:3,1:3]
target.exps <- tcga_expr[symbol,]
other.expr <- tcga_expr[-which(rownames(tcga_expr)==symbol),]
sgcor <- as.data.frame(cor(t(other.expr), t(target.exps)))
colnames(sgcor) <- "r_pearson"
sgcor$pval_pearson <- apply(other.expr, 1, function(x) (cor.test(x, t(target.exps))$p.value))
sgcor_spearman <- as.data.frame(cor(t(other.expr), t(target.exps), method = "spearman"))
colnames(sgcor_spearman) <- "r_spearman"
sgcor_spearman$pval_spearman <- apply(other.expr, 1, function(x)(cor.test(x, t(target.exps), method = "spearman")$p.value))
cors <- cbind(sgcor, sgcor_spearman)
cors$gene <- rownames(other.expr)
head(cors)
write.table(cors,file="1.allCOR.txt",sep="\t",col.names = NA, row.names=T,quote=F) 
newcor <- cors[!(is.na(cors$r_pearson)),]
dim(newcor)
sortcor <- newcor[order(newcor$r_pearson, newcor$r_spearman, decreasing = T),]
topcor <- sortcor[c(1:topnumber, 
                    (nrow(sortcor) - topnumber + 1):nrow(sortcor)),] 
rownames(topcor)
genes <- c(symbol,rownames(topcor))
genesexps <- as.data.frame(t(tcga_expr[genes,]))
sortgenesexps <- genesexps[order(genesexps[,1]),]
write.table(sortgenesexps,file="2.topGeneExp.txt",sep="\t",col.names = NA, row.names=T,quote=F) 
samplenum <- nrow(sortgenesexps)
if(is.null(pointnumber)){
  pointnumber=samplenum
}
group <- as.integer(cut(1:samplenum, breaks=c(seq(from=0.5, to=samplenum, by=samplenum/pointnumber), samplenum+0.5)))
ddf <- data.frame(row.names = 1:pointnumber)
for( i in 1:(1 + topnumber*2)){
  ddf <- cbind(ddf,tapply(sortgenesexps[,i],group,median))
}
colnames(ddf) <- c(symbol,topcor$gene)
mddf <- melt(ddf,id.vars=symbol)
mddf$r <- topcor[mddf$variable,]$r_spearman
mddf$P <- topcor[mddf$variable,]$pval_spearman 
plist <- dlply(mddf, .(variable), function(trig){ggplot(trig, aes_string(x=symbol, y="value")) +
    geom_point() +
    ylab(unique(trig$variable)) +
    theme_classic() +
    ggtitle(paste0("r = ", round(unique(trig$r),2),
                   "\nP = ", sprintf("%1.1e", unique(trig$P)))) +
    geom_smooth(method = "lm", se=F, colour = "#206BB5")})

pg <- plot_grid(plotlist = plist, ncol=1, align = "hv")
pg
filename <- paste0(symbol, "_cor.pdf")
ggsave(filename, width = 4, height = 4)