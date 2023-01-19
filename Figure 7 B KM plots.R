

#setwd("C:\\Users\\M\\Desktop\\um\\23.geneSurvival")         #工作目录（需修改）
library(survminer)
library(survival)
rt=read.table("logforKM.txt",header=T,sep="\t",check.names=F,row.names=1)      #读取文件
#rt$futime=rt$futime/365                                               #如果以月为单位，除以30；以年为单位，除以365


#rt[,4:ncol(rt)] <- as.data.frame(scale(rt[,4:ncol(rt)]))

#rt[,4:ncol(rt)] <- as.data.frame(log2(rt[,4:ncol(rt)]+1.7))
#rt[,4:ncol(rt)] <- as.data.frame(scale(rt[,4:ncol(rt)]))

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
                         #risk.table=TRUE,
                         risk.table.col = "strata", # Change risk table color by groups, 
                         #linetype = "strata", # Change line type by groups
                         #surv.median.line = "hv", # Specify median survival
                         
                         pval =TRUE,
                         conf.int =F,
                         xlab ='Time in years', 
                         legend.title = "Level",
                         legend.labs = c("High", "Low"),
                         
                         #ggtheme = theme_light(), #ggtheme = theme_bw() ggtheme = theme_minimal() ggtheme = theme_light()
                         
    )
    
    print(survp, newpage = FALSE)
    dev.off()
  }
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)

