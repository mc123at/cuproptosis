

#setwd("C:\\Users\\M\\Desktop\\um\\23.geneSurvival")         #工作目录（需修改）
library("survminer")
library(survival)
rt=read.table("sel.immune.LUAD.txt",header=T,sep="\t",check.names=F,row.names=1)      #读取文件
#rt$futime=rt$futime/365                                               #如果以月为单位，除以30；以年为单位，除以365
outTab=data.frame()

for(gene in colnames(rt[,3:ncol(rt)])){
  # a=rt[,gene]
  
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
    if(pValue<0.001){
      pValue="<0.001"
    }else{
      pValue=round(pValue,3)
      pValue=paste0("=",pValue)
    }
    pdf(file=paste(gene,".survival.pdf",sep=""), height = 2.5, width =2)
    survp  <- ggsurvplot(fit,palette = c('#FF69B4', '#87CEFA'), 
                         title = paste("", gene, sep="") , 
                         
                         censor=F,
                         #  risk.table=TRUE,
                         risk.table.col = "strata", # Change risk table color by groups, 
                         #linetype = "strata", # Change line type by groups
                         # surv.median.line = "hv", # Specify median survival
                         break.time.by = 5,
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

