library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor



bioRISK=function(inputFile=null,outFile=null){
  
  
  data=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)  
  #data <- read.csv("easy_input_risk.csv", row.names = 1)
  data[1:2, 1:4]
  
  bestvars <- read.table("geneCoef.txt")$V1
  bestvars
  
  # risk score，用于画顶部散点图
  rs <- data$riskScore
  names(rs) <- rownames(data)
  rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
  # 用中值分组
  rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
  head(rs_data)
  
  # follow-up，用于画中间B图
  surv_data <- data.frame(x=1:length(rs),
                          t=data[names(sort(rs)),'futime']*12,
                          s=data[names(sort(rs)),'fustat']) 
  surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
  head(surv_data)
  
  # 提取signature对应的data，并按risk score排序，用于画底部热图
  exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]
  exp_data[1:2,1:4]
  
  #开始画图
  #A - risk score
  plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
    geom_point(aes(col=Risk),size=1)+
    scale_color_manual(labels=c("High-risk","Low-risk"), 
                       #guide_legend(guide = NULL), #如果不想画图例就删掉#
                       name="Risk score", values =c("#FF69B4", "#87CEFA")) + 
    
    # 画竖向虚线
    geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                     y = 0, 
                     xend = sum(rs_data$Risk=="Low-risk"), 
                     yend = max(rs_data$rs)), linetype="dashed", size = 1)+
    # 画横线
    #geom_segment(aes(x=0,y=median(rs_data$rs),
    #                 xend=nrow(rs_data),
    #                 yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
    
    # 写文字Cutoff:
    #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
    #            y=median(rs_data$rs)+8,
    #            label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
    #        col ="black",size = 4,alpha=0.8)+
    
  theme(axis.title.x=element_blank()) +
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
    labs(y="Risk score",x="",fill="Risk") +
    #scale_colour_discrete(name="Risk scores") +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
          axis.text.x=element_blank())
  
  #plot.A
  
  #B - follow-up
  plot.B <- ggplot(surv_data,aes(x=x,y=t))+
    geom_point(aes(col=Status),shape=2,size=0.3)+
    geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=1,linetype="dashed")+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
    scale_color_manual(labels=c("Alive","Dead"),
                       values =c("#87CEFA","#FF69B4"))+
    labs(y="OS (months)",x="")+
    theme_classic()+
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
          axis.text.x=element_blank())
  #plot.B
  
  #C - signature
  tmp <- t(scale(exp_data))
  tmp[tmp > 2] = 2
  tmp[tmp < -2] = -2
  reorder_cormat <- function(cormat){
    dd <- dist(cormat)
    hc <- hclust(dd,method = "average")
    cormat <-cormat[hc$order,]
  }
  tmp1 <- reorder_cormat(tmp)
  tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low-risk",-2,2))
  tmp.m <- melt(tmp1)
  
  p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
    geom_tile(aes(fill = value)) 
  
  plot.C <- p2 + scale_fill_gradient2(name="Exp", low="#87CEFA", high="#FF69B4", mid="white") +
    labs(x = "", y = "")+
    theme_classic()+
    theme(legend.title = element_text(size = 12), legend.position = "right",
          axis.line = element_blank(),
          axis.ticks=element_blank(),
          axis.text.x=element_blank())
  #plot.C
  
  
  #拼图
  plot_grid(plot.A, plot.B, plot.C,
            labels = c("", "",""), # 或者按顺序标注ABC
            rel_heights = c(0.5,1,1), # 3个图的比例
            #label_x=0,
            #label_y=1,
            align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)
  
  ggsave(outFile, width = 3.5, height = 5.8)
}



bioRISK(inputFile="RISK.1.LUAD_Risk.txt",outFile="1.pdf")
bioRISK(inputFile="RISK.2.GSE72094.txt",outFile="2.pdf")
bioRISK(inputFile="RISK.3.GSE68465.txt",outFile="3.pdf")

#bioRISK(inputFile="RISK.3.GSE68465.txt",outFile="3.pdf")
#bioRISK(inputFile="RISK.4.train.expPFI.5years.txt",outFile="4.pdf")


#bioRISK(inputFile="RISK.5.GSE72094.ExpOs.txt",outFile="5.pdf")
#bioRISK(inputFile="RISK.6.GSE72094.ExpOs.5years.txt",outFile="6.pdf")

#bioRISK(inputFile="RISK.7.GSE68465ExpTime.txt",outFile="7.pdf")
#bioRISK(inputFile="RISK.8.GSE68465.OS.5years.txt",outFile="8.pdf")