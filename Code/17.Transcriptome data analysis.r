library(tidyr)
library(dplyr)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

setwd('/disk213/xieqq/JINHUA138/Transcriptome_analysis.3.RNAseq_new')
Sample <- read.csv("gene_log-tpm-mean_matrix.csv", row.names=1)
#monocle <- read.csv("/disk213/xieqq/sc/Monocle3/Trajectory_genes.csv")
monocle <- read.csv("/disk213/xieqq/JINHUA138.sc/Lamian/trajectory_tests/Trajectory_Genes.csv")
gene=intersect(rownames(Sample), monocle$gene)
Sample <- subset(Sample, rownames(Sample) %in% gene)
monocle <- subset(monocle, monocle$gene %in% intersect(rownames(Sample), monocle$gene))
JH <- Sample[,grepl("JH", colnames(Sample))]
LY <- Sample[,grepl("LY", colnames(Sample))]

#step 1
setwd('./time_warp')
data <- data.frame(gene=monocle$gene,segment=monocle$segment)
data$segment[which(data$segment=="duodenum")]="DU"
data$segment[which(data$segment=="jejunum")]="JE"
data$segment[which(data$segment=="ileum")]="IL"
data$segment[which(data$segment=="cecum")]="CE"
data$segment[which(data$segment=="colon")]="CO"
data <- data[!duplicated(data), ]

judgeJH <- NULL
judgeLY <- NULL
for (i in c(1:nrow(data))) {
  gene <- data$gene[i]
  segment <- data$segment[i]
  
  newJH <- JH[which(rownames(JH)==gene),grepl(segment, colnames(JH))]
  testJH <- data.frame(gene=gene,segment=segment,title=paste(gene,segment,sep="_"),T0_60=NA,T60_90=NA,T90_180=NA,T180_240=NA)
  if (newJH[,1]>newJH[,2]){testJH$T0_60="1"}else{testJH$T0_60="0"}
  if (newJH[,2]>newJH[,3]){testJH$T60_90="1"}else{testJH$T60_90="0"}
  if (newJH[,3]>newJH[,4]){testJH$T90_180="1"}else{testJH$T90_180="0"}
  if (newJH[,4]>newJH[,5]){testJH$T180_240="1"}else{testJH$T180_240="0"}
  judgeJH <- rbind(judgeJH,testJH)
  
  newLY <- LY[which(rownames(LY)==gene),grepl(segment, colnames(LY))]
  testLY <- data.frame(gene=gene,segment=segment,title=paste(gene,segment,sep="_"),T0_60=NA,T60_90=NA,T90_180=NA) 
  if (newLY[,1]>newLY[,2]){testLY$T0_60="1"}else{testLY$T0_60="0"}
  if (newLY[,2]>newLY[,3]){testLY$T60_90="1"}else{testLY$T60_90="0"}
  if (newLY[,3]>newLY[,4]){testLY$T90_180="1"}else{testLY$T90_180="0"}
  judgeLY <- rbind(judgeLY,testLY)
}
judgeJH$JH <- paste(judgeJH$T0_60,judgeJH$T60_90,judgeJH$T90_180,judgeJH$T180_240,sep="-")
judgeLY$LY <- paste(judgeLY$T0_60,judgeLY$T60_90,judgeLY$T90_180,sep="-")
judge <- inner_join(judgeJH[,c(1,2,3,8)],judgeLY[,c(1,2,3,7)],by=c("gene","segment","title"))

#step 2: match,warp(forward/backward),out
match_list <- list()
match_list[[1]] <- c("0-0-0-0","0-0-0")
match_list[[2]] <- c("0-0-0-1","0-0-0")
match_list[[3]] <- c("0-0-1-0","0-0-1")
match_list[[4]] <- c("0-0-1-1","0-0-1")
match_list[[5]] <- c("0-1-0-0","0-1-0")
match_list[[6]] <- c("0-1-0-1","0-1-0")
match_list[[7]] <- c("0-1-1-0","0-1-1")
match_list[[8]] <- c("0-1-1-1","0-1-1")
match_list[[9]] <- c("1-0-0-0","1-0-0")
match_list[[10]] <- c("1-0-0-1","1-0-0")
match_list[[11]] <- c("1-0-1-0","1-0-1")
match_list[[12]] <- c("1-0-1-1","1-0-1")
match_list[[13]] <- c("1-1-0-0","1-1-0")
match_list[[14]] <- c("1-1-0-1","1-1-0")
match_list[[15]] <- c("1-1-1-0","1-1-1")
match_list[[16]] <- c("1-1-1-1","1-1-1")

judge_match <- NULL
for (i in c(1:16)){
  out <- judge[which(judge$JH==match_list[[i]][1]&judge$LY==match_list[[i]][2]),]
  judge_match <- rbind(judge_match,out)
}
count_match <- nrow(judge_match)
genecount_match <- length(unique(judge_match$gene))

forward_list <- list()
forward_list[[1]] <- c("0-0-0-1","0-0-1")
forward_list[[2]] <- c("0-0-0-1","0-1-1")
forward_list[[3]] <- c("0-0-1-0","0-1-0")
forward_list[[4]] <- c("0-0-1-1","0-1-1")
forward_list[[5]] <- c("0-1-1-0","0-1-0")
forward_list[[6]] <- c("1-0-0-1","1-0-1")
forward_list[[7]] <- c("1-1-0-0","1-0-0")
forward_list[[8]] <- c("1-1-0-1","1-0-1")
forward_list[[9]] <- c("1-1-1-0","1-1-0")
forward_list[[10]] <- c("1-1-1-0","1-0-0")

judge_forward <- NULL
for (i in c(1:10)){
  out <- judge[which(judge$JH==forward_list[[i]][1]&judge$LY==forward_list[[i]][2]),]
  judge_forward <- rbind(judge_forward,out)
}
count_forward <- nrow(judge_forward)
genecount_forward <- length(unique(judge_forward$gene))

backward_list <- list()
backward_list[[1]] <- c("0-0-1-1","0-0-0")
backward_list[[2]] <- c("0-1-1-1","0-0-1")
backward_list[[3]] <- c("1-0-0-0","1-1-0")
backward_list[[4]] <- c("1-1-0-0","1-1-1")

judge_backward <- NULL
for (i in c(1:4)){
  out <- judge[which(judge$JH==backward_list[[i]][1]&judge$LY==backward_list[[i]][2]),]
  judge_backward <- rbind(judge_backward,out)
}
count_backward <- nrow(judge_backward)
genecount_backward <- length(unique(judge_backward$gene))

judge_warp <- rbind(judge_forward,judge_backward)
count_warp <- nrow(judge_warp)
genecount_warp <- length(unique(judge_warp$gene))

filter <- union(judge_forward$title, judge_backward$title)
filter <- union(filter, judge_match$title)
judge_out <- judge[!(judge$title %in% filter), ]
count_out <- nrow(judge_out)
genecount_out <- length(unique(judge_out$gene))

judge_forward$type <- "forward"
judge_backward$type <- "backward"
judge_match$type <- "match"
judge_out$type <- "out"
judge <- rbind(judge_forward[,-3],judge_backward[,-3],judge_match[,-3],judge_out[,-3])
write.table(judge, "cluster.txt", row.names=F, quote=F) 

count <- matrix(data=NA, nrow=2, ncol=5, dimnames=list(c("count","genecount"),c("forward","backward","warp","match","out")))
count[1,] <- c(count_forward,count_backward,count_warp,count_match,count_out)
count[2,] <- c(genecount_forward,genecount_backward,genecount_warp,genecount_match,genecount_out)
write.table(count, "cluster_count.txt", row.names=T, quote=F)

#step 3
library(grid)
library(reshape2)
count <- read.table("cluster_count.txt",header=T)
count$group <- rownames(count)
count <- melt(count)
count$variable <- factor(count$variable,levels=c("out","match","warp","forward","backward"))
count1 <- subset(count,variable %in% c("warp","match","out"))
count2 <- subset(count,variable %in% c("forward","backward"))

p1 <- ggplot(count1,aes(x=variable, y=value, fill=group, label=value))+
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.6)+
  geom_text(size=6,vjust=-0.5,position = position_dodge(0.7))+
  scale_fill_manual(values=c("#D95B5B","#9191C8"))+
  labs(x="", y="Count", title="")+
  theme_bw()+ 
  theme(axis.text.x=element_text(color="black",size=16),
        axis.text.y=element_text(color="black",size=16),
        axis.title.x=element_text(color="black",size=18),
        axis.title.y=element_text(color="black",size=18),
        legend.text=element_text(color="black",size=16),
        legend.title=element_text(color="black",size=18),
        legend.position=c(0.5, 0.95), 
        legend.direction = "horizontal", 
        axis.line = element_line(colour = "black"), 
        panel.border = element_blank(), 
        panel.grid.major=element_blank(),   
        panel.grid.minor=element_blank())

p2 <- ggplot(count2,aes(x=variable, y=value, fill=group, label=value))+
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.6)+
  geom_text(size=4,vjust=-0.5,position = position_dodge(0.7))+
  scale_fill_manual(values=c("#D95B5B","#9191C8"))+
  labs(x="", y="", title="")+
  theme_bw()+ 
  theme(axis.text.x=element_text(color="black",size=16),
        axis.text.y=element_text(color="black",size=16),
        axis.title.x=element_text(color="black",size=18),
        axis.title.y=element_text(color="black",size=18),
        legend.text=element_text(color="black",size=16),
        legend.title=element_text(color="black",size=18),
        legend.position="none",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
pdf(file=paste("cluster_count.pdf"), width=8, height=8)
print(p1)
vie <- viewport(width=0.45,height=0.4,x=0.75,y=0.7)
print(p2,vp=vie)
dev.off()

#step 4
library(ggpubr)
library(cowplot)
inf_list <- list()
inf_list[[1]] <- c("DU","#E76F51")
inf_list[[2]] <- c('JE','#F4A261')
inf_list[[3]] <- c('IL','#E9C46A')
inf_list[[4]] <- c('CE','#2A9D8F')
inf_list[[5]] <- c('CO','#264653')

my_plot <- function(data,color,ylab) {
  if (min(data$value)%%1>0.5){y_min <- floor(min(data$value))+0.5} else{y_min <- floor(min(data$value))}
  if (max(data$value)%%1<0.5){y_max <- ceiling(max(data$value))-0.5} else{y_max <- ceiling(max(data$value))}
  p1 <- ggplot(data,aes(x=variable, y=value, group=gene))+
    geom_line(color="gray90",linewidth=0.8) + 
    geom_hline(yintercept=0,linetype=2) +
    stat_summary(aes(group=1), fun=mean, geom="line", linewidth=1.2, color=color) + 
    labs(x="", y=ylab, title="")+
    scale_y_continuous(limits=c(y_min, y_max))+
    theme_bw()+ 
    theme(axis.text.x=element_text(color="black",size=14),
          axis.text.y=element_text(color="black",size=14),
          axis.title.x=element_text(color="black",size=16),
          axis.title.y=element_text(color="black",size=16),
          legend.text=element_text(color="black",size=14),
          legend.title=element_text(color="black",size=16),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank())
  return(p1)}

cluster_all <- read.table("cluster.txt",header=T)
cluster <- cluster_all[which(cluster_all$type=="forward"),]
cluster$group <- paste(cluster$JH,cluster$LY,sep="-")
cluster$cluster <- NA
for (i in c(1:length(unique(cluster$group)))){
  cluster$cluster[which(cluster$group==unique(cluster$group)[i])]=paste("Cluster",i,sep="-")
}
new_cluster <- left_join(cluster_all[,c(1,2,5)],cluster[,c(1,2,5,7)],by=c('gene', 'segment', 'type'))
new_cluster$type <- factor(new_cluster$type,levels=c("forward","backward","match","out"),labels=c("warp (forward)","warp (backward)","match","out"))
new_cluster$segment <- factor(new_cluster$segment,levels=c("DU","JE","IL","CE","CO"),labels=c("Duodenum","Jejunum","Ileum","Cecum","Colon"))
new_cluster$cluster <- factor(new_cluster$cluster,levels=c("Cluster-1","Cluster-2","Cluster-3","Cluster-4","Cluster-5","Cluster-6","Cluster-7","Cluster-8","Cluster-9","Cluster-10"))
new_cluster <- arrange(new_cluster,type,segment,cluster,gene)
write.csv(new_cluster, "cluster_all.csv", row.names=F, quote=F)

for (s in c(1:5)){
  segment=inf_list[[s]][1]
  file_1 <- JH[,grepl(segment, colnames(JH))]
  file_2 <- LY[,grepl(segment, colnames(LY))]
  cluster_new <- cluster[which(cluster$segment==segment),]
  for (i in unique(cluster_new$cluster)){
    gene <- cluster_new$gene[which(cluster_new$cluster==i)]
    file_JH <- subset(file_1, rownames(file_1) %in% gene)
    colnames(file_JH) <- c("0d","60d","90d","180d","240d")
    file_JH$gene <- rownames(file_JH)
    file_JH = melt(file_JH)
    ylab <- "Normalized mean value"
    file_LY <- subset(file_2, rownames(file_2) %in% gene)
    colnames(file_LY) <- c("0d","60d","90d","180d")
    file_LY$gene <- rownames(file_LY)
    file_LY = melt(file_LY)
    text = paste(segment,i,sep="-")
    P1 <- my_plot(file_JH, inf_list[[s]][2],ylab)
    P2 <- my_plot(file_LY, inf_list[[s]][2],"")
    P3 <- ggdraw()+
      draw_plot(P1,x=0,y=0,width=0.7,height=0.7)+
      draw_plot(P2,x=0.5,y=0.5,width=0.5,height=0.5)+
      draw_plot_label(label=c("JH","LY"),fontface="plain",size=16,x=c(0.1,0.6),y=c(0.7,1))+
      draw_plot_label(label=c(text),fontface="plain",color=inf_list[[s]][2],size=16,x=c(0.55,0.85),y=c(0.15,0.65))
    pdf(file=paste0(segment,"_",i,".pdf"), width=10, height=10)
    print(P3)
    #print(ggarrange(P1,P2, nrow=2, common.legend=TRUE, legend="right"))
    dev.off()
  }
}
