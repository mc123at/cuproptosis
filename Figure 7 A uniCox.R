library(survival) 
library(glmnet) 
library(forestplot) 
library(randomForestSRC)
library(dplyr)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
display.progress = function (index, totalN, breakN=20) {
        if ( index %% ceiling(totalN/breakN)  ==0  ) {
                cat(paste(round(index*100/totalN), "% ", sep=""))
        }
}  
tcga.surv <- read.table("CIBERSORT-Results.downFromGDC for plot.txt",sep="\t",row.names = 1, check.names = F, stringsAsFactors = F, header = T)
head(tcga.surv)
tcga.surv[,3:ncol(tcga.surv)] <- as.data.frame(scale(tcga.surv[,3:ncol(tcga.surv)]))
write.table(tcga.surv,file="logforKM.txt",sep="\t",row.names=T,col.names=NA,quote=F)
uni.tcga=data.frame()
for(i in colnames(tcga.surv[,3:ncol(tcga.surv)])){
  cox <- coxph(Surv(futime, fustat) ~ tcga.surv[,i], data = tcga.surv)
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
write.table(uni.tcga,file="UnicoxFororder.txt",sep="\t",row.names=F,quote=F)
uni.tcga <- read.table("UnicoxFororder.txt",sep="\t",check.names = F, stringsAsFactors = F, header = T)
hrtable <- rbind(uni.tcga)
tabletext <- cbind(c("Variable",hrtable$variable),
                   c("HRHRHRHR",format(round(as.numeric(hrtable$HR),3),nsmall = 3)),
                   c("Lower95%CI",format(round(as.numeric(hrtable$lower.95CI),3),nsmall = 3)),
                   c("Upper95%CI",format(round(as.numeric(hrtable$upper.95CI),3),nsmall = 3)),
                   c("PPPPPP-value",formatC(as.numeric(hrtable$p), format = "e", digits = 2)))
tabletext
nrow(tabletext) + 1 
pdf("forestplot of risk table.pdf", width = 12, height = 10)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))), 
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(-0.5,0,0.5),
           lwd.xaxis=2,
           xlab=expression("log"[2]~"HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="black"),
                           "24" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(.75,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())