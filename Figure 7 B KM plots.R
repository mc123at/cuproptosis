library(survminer)
library(survival)
rt=read.table("logforKM.txt",header=T,sep="\t",check.names=F,row.names=1)
outTab=data.frame()
for(gene in colnames(rt[,3:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)
  if(pValue<1.05){
    if(pValue<0.001){
      pValue="<0.001"
    }else{
      pValue=round(pValue,3)
      pValue=paste0("=",pValue)
    }
    pdf(file=paste(gene,".survival.pdf",sep=""), height = 2.5, width = 2.5)
    survp  <- ggsurvplot(fit,palette = c('#FF69B4', '#87CEFA'), 
                         title = paste("", gene, sep="") , 
                         censor=T,
                         break.time.by = 5,
                         risk.table.col = "strata", 
                         pval =TRUE,
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