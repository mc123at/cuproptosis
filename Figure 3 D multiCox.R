library(survival) 
library(glmnet) 
library(forestplot) 
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
display.progress = function (index, totalN, breakN=20) {
        if ( index %% ceiling(totalN/breakN)  ==0  ) {
                cat(paste(round(index*100/totalN), "% ", sep=""))
        }
}  
tcga.surv <- read.table("1.LUADcox.clinic.riskscore.txt", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
#head(tcga.surv)
mulcox.tcga <- summary(coxph(Surv(OS.time, OS) ~ ., data = tcga.surv))
mulcox.tcga <- data.frame(variable = rownames(mulcox.tcga$conf.int),
                          HR = mulcox.tcga$conf.int[,1],
                          lower.95CI = mulcox.tcga$conf.int[,3],
                          upper.95CI = mulcox.tcga$conf.int[,4],
                          p = mulcox.tcga$coefficients[,5],
                          stringsAsFactors = F)
rownames(mulcox.tcga) <- NULL
geo.surv <- read.table("2.GSE72094cox.clinic.riskscore.txt", row.names = 1, check.names = F, stringsAsFactors = F, header = T)

mulcox.geo <- summary(coxph(Surv(OS.time, OS) ~ ., data = geo.surv))
mulcox.geo <- data.frame(variable = rownames(mulcox.geo$conf.int),
                         HR = mulcox.geo$conf.int[,1],
                         lower.95CI = mulcox.geo$conf.int[,3],
                         upper.95CI = mulcox.geo$conf.int[,4],
                         p = mulcox.geo$coefficients[,5],
                         stringsAsFactors = F)
rownames(mulcox.geo) <- NULL
geo.surv22 <- read.table("3.GSE68465cox.clinic.riskscore.txt", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
mulcox.geo22 <- summary(coxph(Surv(OS.time, OS) ~ ., data = geo.surv22))
mulcox.geo22 <- data.frame(variable = rownames(mulcox.geo22$conf.int),
                         HR = mulcox.geo22$conf.int[,1],
                         lower.95CI = mulcox.geo22$conf.int[,3],
                         upper.95CI = mulcox.geo22$conf.int[,4],
                         p = mulcox.geo22$coefficients[,5],
                         stringsAsFactors = F)
rownames(mulcox.geo22) <- NULL
hrtable <- rbind(c("TCGA-LUAD",NA,NA,NA,NA),
                 mulcox.tcga,
                 c("GSE72094",NA,NA,NA,NA),
                 mulcox.geo,
                 c("GSE68465",NA,NA,NA,NA),
                 mulcox.geo22)
tabletext <- cbind(c("Variable",hrtable$variable),
                   c("HRHRHRHR",format(round(as.numeric(hrtable$HR),3),nsmall = 3)),
                   c("Lower95%CI",format(round(as.numeric(hrtable$lower.95CI),3),nsmall = 3)),
                   c("Upper95%CI",format(round(as.numeric(hrtable$upper.95CI),3),nsmall = 3)),
                   c("PPPPPP-value",formatC(as.numeric(hrtable$p), format = "e", digits = 2)))
tabletext
nrow(tabletext) + 1 
tabletext[2,] <- c("TCGA-LUAD",NA,NA,NA,NA) 
tabletext[14,] <- c("GSE72094",NA,NA,NA,NA) 
tabletext[28,] <- c("GSE68465",NA,NA,NA,NA) 
pdf("multi.forestplot of risk table.pdf", width = 15, height = 15)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))), 
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#ffc403", lines="#e0ac00", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(-2,0,2,4),
           lwd.xaxis=2,
           xlab=expression("log"[2]~"HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "14" = gpar(lwd=1, col="grey50", lty=2),
                           "28" = gpar(lwd=1, col="grey50", lty=2),
                           "35" = gpar(lwd=2, col="black")),
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