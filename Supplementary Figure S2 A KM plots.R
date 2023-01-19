library("survminer")
library(survival)
rt=read.table("RISK.1.LUAD_Risk.txt",header=T,sep="\t",check.names=F,row.names=1)
outTab=data.frame()
for(gene in colnames(rt[,3:ncol(rt)])){
  res.cut <- surv_cutpoint(rt, time = "futime", event = "fustat",
                           variables = gene)
  res.cat <- surv_categorize(res.cut)
a=res.cat[,gene]
   diff=survdiff(Surv(futime, fustat) ~ a,data = res.cat)
   pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  fit <- survfit(Surv(futime, fustat) ~ a, data = res.cat)
  summary(fit)
  if(pValue<1.05){
    if(pValue<0.000000001){
      pValue="<0.001"
    }else{
      pValue=signif(pValue,3)
      pValue=format(pValue, scientific = TRUE)
    }
    pdf(file=paste(gene,".survival.pdf",sep=""), height = 2.5, width =2)
    survp  <- ggsurvplot(fit,palette = c('#FF69B4', '#87CEFA'), 
                         title = paste("", gene, sep="") , 
                         censor=F,
                         risk.table.col = "strata", 
                         break.time.by = 5,
                         pval =paste0("p = ", pValue),
                         conf.int =F,
                         xlab ='Time in years', 
                         legend.title = "Level",
                         legend.labs = c("High", "Low"),
    )
    
    print(survp, newpage = FALSE)
    dev.off()
  }
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)