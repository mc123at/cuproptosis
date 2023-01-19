#install.packages("timeROC")

library(timeROC)
library(survival)

Sys.setenv(LANGUAGE = "en") #??ʾӢ?ı?????Ϣ
options(stringsAsFactors = FALSE) #??ֹchrת??factor

tr=read.table("1.LUADcox.clinic.riskscore.txt",header=T,sep="\t",check.names=F,row.names=1) 
head(tr)
dim(tr)
ROC.riskscore1 <- timeROC(T=tr$futime, 
                         delta=tr$fustat, marker=tr$RiskScore,
                         #other_markers=as.matrix(tr[,c("age","sex")]),
                         cause=1,
                         weighting="marginal",
                         times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                         #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                         # iid=TRUE,
                         )






ROC.riskscore <- timeROC(T=tr$futime, 
                 delta=tr$fustat, marker=tr$RiskScore,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                 #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                 # iid=TRUE
                 )

ROC.age <- timeROC(T=tr$futime, 
                     delta=tr$fustat, marker=tr$Age,
                     #other_markers=as.matrix(tr[,c("age","sex")]),
                     cause=1,
                     weighting="marginal",
                   times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                   #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                     # iid=TRUE
                     )

ROC.gender <- timeROC(T=tr$futime, 
                         delta=tr$fustat, marker=tr$Gender,
                         #other_markers=as.matrix(tr[,c("age","sex")]),
                         cause=1,
                         weighting="marginal",
                      times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                      #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                         #iid=TRUE,
                     )


ROC.Race <- timeROC(T=tr$futime, 
                          delta=tr$fustat, marker=tr$Race,
                          #other_markers=as.matrix(tr[,c("age","sex")]),
                          cause=1,
                          weighting="marginal",
                    times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                    #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                    #iid=TRUE
                    )

ROC.ethnicity <- timeROC(T=tr$futime, 
                          delta=tr$fustat, marker=tr$Ethnicity,
                          #other_markers=as.matrix(tr[,c("age","sex")]),
                          cause=1,
                          weighting="marginal",
                         times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                         #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                          #iid=TRUE
)


ROC.TumorStage <- timeROC(T=tr$futime, 
                         delta=tr$fustat, marker=tr$TumorStage,
                         #other_markers=as.matrix(tr[,c("age","sex")]),
                         cause=1,
                         weighting="marginal",
                         times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                         #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                         #iid=TRUE
                         )




ROC.PriorMalignancy <- timeROC(T=tr$futime, 
                         delta=tr$fustat, marker=tr$PriorMalignancy,
                         #other_markers=as.matrix(tr[,c("age","sex")]),
                         cause=1,
                         weighting="marginal",
                         times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                         #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                        #iid=TRUE
                        )





ROC.TissueOrOrganOfOrigin <- timeROC(T=tr$futime, 
                         delta=tr$fustat, marker=tr$TissueOrOrganOfOrigin,
                         #other_markers=as.matrix(tr[,c("age","sex")]),
                         cause=1,
                         weighting="marginal",
                         times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                         #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                         #iid=TRUE
                         )


ROC.TobaccoSmokingHistory <- timeROC(T=tr$futime, 
                                     delta=tr$fustat, marker=tr$TobaccoSmokingHistory,
                                     #other_markers=as.matrix(tr[,c("age","sex")]),
                                     cause=1,
                                     weighting="marginal",
                                     times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                                     #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                                     #iid=TRUE
)

ROC.RiskScorePlusStage <- timeROC(T=tr$futime, 
                                     delta=tr$fustat, marker=tr$RiskScorePlusStage,
                                     #other_markers=as.matrix(tr[,c("age","sex")]),
                                     cause=1,
                                     weighting="marginal",
                                     times=c(1,1.5,2,2.5,3,3.5,4,4.5,5),
                                     #times=c(1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                                     #iid=TRUE
)



pdf("2.timeROC.pdf", 5, 4)
plotAUCcurve(ROC.riskscore, conf.int=FALSE, col="red")
plotAUCcurve(ROC.age, conf.int=FALSE, col="darkblue", add=TRUE)
plotAUCcurve(ROC.gender, conf.int=FALSE, col="darkgreen", add=TRUE)
plotAUCcurve(ROC.Race, conf.int=FALSE, col="yellow", add=TRUE)
plotAUCcurve(ROC.ethnicity, conf.int=FALSE, col="pink", add=TRUE)
plotAUCcurve(ROC.TumorStage, conf.int=FALSE, col="green", add=TRUE)
plotAUCcurve(ROC.PriorMalignancy, conf.int=FALSE, col="#9B870C", add=TRUE)


plotAUCcurve(ROC.TissueOrOrganOfOrigin, conf.int=FALSE, col="#577aa1", add=TRUE)
plotAUCcurve(ROC.TobaccoSmokingHistory, conf.int=FALSE, col="#5c1228", add=TRUE)
plotAUCcurve(ROC.RiskScorePlusStage, conf.int=FALSE, col="#0022ff", add=TRUE)
#plotAUCcurve(ROC.riskall, conf.int=FALSE, col="gold", add=TRUE)

# ͼ??????
legend("topright", c("risk score","age","gender","race","ethnicity","TumorStage","PriorMalignancy","TissueOrOrganOfOrigin","TobaccoSmokingHistory","RiskScorePlusStage"),
       col=c("red","darkblue","darkgreen","yellow","pink","green","#9B870C","#577aa1","#5c1228","#0022ff"),
        lty=1, lwd=2, cex=0.8)
dev.off()


#??roc 1,3,5years
pdf("1.ROC.pdf", 5, 4)
plot(ROC.riskscore1,time=1)        
plot(ROC.riskscore1,time=3,add=TRUE,col="blue") 
plot(ROC.riskscore1,time=5,add=TRUE,col="gold") 
legend("bottomright",c("1 year(AUC=      )","3 years(AUC=      )","5 years(AUC=      )"),col=c("red","blue","gold"),lty=1,lwd=2)
dev.off()

ROC.riskscore1
ROC.riskscore1[["AUC"]]

ROC.TumorStage[["AUC"]]

ROC.RiskScorePlusStage[["AUC"]]