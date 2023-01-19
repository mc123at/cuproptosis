library(survival) # 生存分析
library(glmnet) # LASSO回归
library(forestplot) # 绘制森林图
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
display.progress = function (index, totalN, breakN=20) {
        if ( index %% ceiling(totalN/breakN)  ==0  ) {
                cat(paste(round(index*100/totalN), "% ", sep=""))
        }
}  
#train.unicox
tcga.surv <- read.table("1.LUADcox.clinic.riskscore.txt", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
head(tcga.surv)


uni.tcga=data.frame()
for(i in colnames(tcga.surv[,3:ncol(tcga.surv)])){
  cox <- coxph(Surv(OS.time, OS) ~ tcga.surv[,i], data = tcga.surv)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  coef<-signif(coxSummary$coefficients[,"coef"], digits=6)
  HR<-signif(coxSummary$conf.int[,"exp(coef)"], digits=6)
  HR95L<-signif(coxSummary$conf.int[,"lower .95"], digits=6)
  HR95H<-signif(coxSummary$conf.int[,"upper .95"], digits=6)
  z<-signif(coxSummary$coefficients[,"z"], digits=6)
  pvalue<-signif(coxSummary$coefficients[,"Pr(>|z|)"], digits=6) 
  
  
  uni.tcga=rbind(uni.tcga,
               cbind(variable=i,
                     HR=paste0(HR),
                     lower.95CI = paste0(HR95L),
                     upper.95CI = paste0(HR95H),
                     p=paste0(pvalue))
  )
}

rownames(uni.tcga) <- NULL
head(uni.tcga)





#valid1.unicox
geo.surv <- read.table("2.GSE72094cox.clinic.riskscore.txt", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
head(geo.surv)

uni.geo1=data.frame()
for(i in colnames(geo.surv[,3:ncol(geo.surv)])){
  cox <- coxph(Surv(OS.time, OS) ~ geo.surv[,i], data = geo.surv)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  coef<-signif(coxSummary$coefficients[,"coef"], digits=6)
  HR<-signif(coxSummary$conf.int[,"exp(coef)"], digits=6)
  HR95L<-signif(coxSummary$conf.int[,"lower .95"], digits=6)
  HR95H<-signif(coxSummary$conf.int[,"upper .95"], digits=6)
  z<-signif(coxSummary$coefficients[,"z"], digits=6)
  pvalue<-signif(coxSummary$coefficients[,"Pr(>|z|)"], digits=6) 
  
  
  uni.geo1=rbind(uni.geo1,
                 cbind(variable=i,
                       HR=paste0(HR),
                       lower.95CI = paste0(HR95L),
                       upper.95CI = paste0(HR95H),
                       p=paste0(pvalue))
  )
}

rownames(uni.geo1) <- NULL
head(uni.geo1)




#valid2.unicox
geo.surv22 <- read.table("3.GSE68465cox.clinic.riskscore.txt", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
head(geo.surv22)

uni.geo22=data.frame()
for(i in colnames(geo.surv22[,3:ncol(geo.surv22)])){
  cox <- coxph(Surv(OS.time, OS) ~ geo.surv22[,i], data = geo.surv22)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  coef<-signif(coxSummary$coefficients[,"coef"], digits=6)
  HR<-signif(coxSummary$conf.int[,"exp(coef)"], digits=6)
  HR95L<-signif(coxSummary$conf.int[,"lower .95"], digits=6)
  HR95H<-signif(coxSummary$conf.int[,"upper .95"], digits=6)
  z<-signif(coxSummary$coefficients[,"z"], digits=6)
  pvalue<-signif(coxSummary$coefficients[,"Pr(>|z|)"], digits=6) 
  
  
  uni.geo22=rbind(uni.geo22,
                 cbind(variable=i,
                       HR=paste0(HR),
                       lower.95CI = paste0(HR95L),
                       upper.95CI = paste0(HR95H),
                       p=paste0(pvalue))
  )
}

rownames(uni.geo22) <- NULL
head(uni.geo22)













#绘制森林图
hrtable <- rbind(c("TCGA-LUAD",NA,NA,NA,NA),
                 uni.tcga,
                 c("GSE72094",NA,NA,NA,NA),
                 uni.geo1,
                 c("GSE68465",NA,NA,NA,NA),
                 uni.geo22)
tabletext <- cbind(c("Variable",hrtable$variable),
                   c("HRHRHRHR",format(round(as.numeric(hrtable$HR),3),nsmall = 3)),
                   c("Lower95%CI",format(round(as.numeric(hrtable$lower.95CI),3),nsmall = 3)),
                   c("Upper95%CI",format(round(as.numeric(hrtable$upper.95CI),3),nsmall = 3)),
                   c("PPPPPP-value",formatC(as.numeric(hrtable$p), format = "e", digits = 2)))
tabletext

nrow(tabletext) + 1 #把这个数字写入hrzl_lines参数的第四行

# 按需设置，因为注意到第二行的NA变成了字符串，因此会显示在最终的森林图里，这里改为NA
tabletext[2,] <- c("TCGA-LUAD",NA,NA,NA,NA) 
# 按需设置，因为注意到第九行的NA变成了字符串，因此会显示在最终的森林图里，这里改为NA
tabletext[14,] <- c("GSE72094",NA,NA,NA,NA) 
tabletext[28,] <- c("GSE68465",NA,NA,NA,NA) 


pdf("uni.forestplot of risk table.pdf", width = 15, height = 15)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),#log2(HR)
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))), #log2(95%置信区间下限)
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),#log2(95%置信区间上限)
           graph.pos=6,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#ffc403", lines="#e0ac00", zero = "black"),#box颜色
           boxsize=0.4,#box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,#不显示区间
           zero=0,#zero线横坐标
           lwd.zero=2,#zero线宽
           xticks = c(-2,0,2,4),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab=expression("log"[2]~"HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第二行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                           "14" = gpar(lwd=1, col="grey50", lty=2),#第九行顶部加灰色虚线
                           "28" = gpar(lwd=1, col="grey50", lty=2),#第九行顶部加灰色虚线
                           "35" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(.75,"cm"),#固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())


########


