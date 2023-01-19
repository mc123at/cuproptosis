

library(survival)
library("survminer")
#setwd("D:\\biowolf\\metabolism\\16.survival")              #设置工作目录

bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")                   #读取输入文件
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     #conf.int=TRUE,
                     pval= paste0("p-value = ", pValue),
                     pval.size=4,
                     risk.table=TRUE,
                     legend.labs=c("high", "low"),
                     legend.title="Risk level",
                     xlab="Time (years)",
                     #break.time.by = 1,
                     surv.median.line="hv",
                     #xlim = c(0, 10),
                     risk.table.title="",
                     palette=c("#FF69B4", "#87CEFA"),
                     ncensor.plot.title = "",
                     ncensor.plot = TRUE,
                     risk.table.height=.25,
  )
  pdf(file=outFile,onefile = FALSE,width = 2.5,height =6)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="RISK.1.LUAD_Risk.txt",outFile="1.pdf")
bioSurvival(inputFile="RISK.2.GSE72094.txt",outFile="2.pdf")
bioSurvival(inputFile="RISK.3.GSE68465.txt",outFile="3.pdf")


