#????????
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
Sys.setenv(LANGUAGE = "en") #??ʾӢ?ı?????Ϣ
options(stringsAsFactors = FALSE) #??ֹchrת??factor

#????????
symbol <- "Riskscore"
topnumber = 10  #ȡ??????????Top X?͸?????????Top X?Ļ???

# plot with special pointnumber by groupping gene expression levels into equally sized bins
pointnumber = 150 
# plot without bins
#pointnumber = NULL

#?????ļ?
#tcga_expr <- read.csv("TMB and GNPNAT1.txt", row.names = 1)

tcga_expr=read.table("selectgenes.Exp.txt", header=T,sep="\t",row.names=1,check.names=F)
tcga_expr[1:3,1:3]

#????????ϵ??
# Ŀ??????
target.exps <- tcga_expr[symbol,]
# ????????
other.expr <- tcga_expr[-which(rownames(tcga_expr)==symbol),]



# pearson
sgcor <- as.data.frame(cor(t(other.expr), t(target.exps))) #Ĭ??pearson
colnames(sgcor) <- "r_pearson"
# ????pvalue??????ʱ???ϳ?
sgcor$pval_pearson <- apply(other.expr, 1, function(x) (cor.test(x, t(target.exps))$p.value))



# spearman
sgcor_spearman <- as.data.frame(cor(t(other.expr), t(target.exps), method = "spearman"))
colnames(sgcor_spearman) <- "r_spearman"
sgcor_spearman$pval_spearman <- apply(other.expr, 1, function(x)(cor.test(x, t(target.exps), method = "spearman")$p.value))






# ??????ϵ????pvalue??д??cors??
cors <- cbind(sgcor, sgcor_spearman)
cors$gene <- rownames(other.expr)
head(cors)
write.table(cors,file="1.allCOR.xls",sep="\t",col.names = NA, row.names=T,quote=F) #??һ??????Ϊ?գ?????NA???Խ???




# ȡ????ϵ??????ǰ??λ?Ļ??????ڲ??????????޸?topnumber
# ?˴?ȡǰX?ͺ?X??Ҳ??????????TopX?͸?????TopX
newcor <- cors[!(is.na(cors$r_pearson)),]
dim(newcor)

sortcor <- newcor[order(newcor$r_pearson, newcor$r_spearman, decreasing = T),]
topcor <- sortcor[c(1:topnumber, #??????Top X
                    (nrow(sortcor) - topnumber + 1):nrow(sortcor)),] #??????Top X
rownames(topcor)

# ??ȡ????ϵ??????ǰ??λ?????ı???????,Ŀ??????Ҳ??д??table??Ϊ??һ?л???
genes <- c(symbol,rownames(topcor))
genesexps <- as.data.frame(t(tcga_expr[genes,]))
sortgenesexps <- genesexps[order(genesexps[,1]),]
write.table(sortgenesexps,file="2.topGeneExp.xls",sep="\t",col.names = NA, row.names=T,quote=F) #??һ??????Ϊ?գ?????NA???Խ???


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

# ??ͼ????ʾpearson??r??P
# ??????ʾspearman??r??P???ͻ???r_spearman??pval_spearman
mddf$r <- topcor[mddf$variable,]$r_pearson
mddf$P <- topcor[mddf$variable,]$pval_pearson 

# mddf$r <- topcor[mddf$variable,]$r_spearman
# mddf$P <- topcor[mddf$variable,]$pval_spearman 










#??ʼ??ͼ

#??��??ͼ

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
ggsave(filename, width =10, height = 9)









