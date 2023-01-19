rt=read.table("CIBERSORT-Results.downFromGDC for plot.txt",sep="\t",header=T,row.names=1,check.names=F)
library(corrplot)
pdf("corHeatmap.pdf",height=13,width=13)  
res1 <- cor.mtest(rt, conf.level = .95)
M = cor(rt)
col3 <- colorRampPalette(c("blue", "white", "red")) 
corrplot(M,
         type = "upper",
         method = "color",
         order = "hclust",
   
         tl.col="black",
         addCoef.col = "black",
         p.mat = res1$p, sig.level = .05,
         
         number.cex = 0.8,
         col=col3(100),
         )
dev.off()