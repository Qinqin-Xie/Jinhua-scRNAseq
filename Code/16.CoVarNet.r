# devtools::install_github(repo = "https://github.com/QiangShiPKU/CoVarNet")
library(CoVarNet)
library(NMF)
library(igraph)
library(dplyr)
library(tibble)
library(grid)
library(ggplot2)
library(tidytext)


setwd("/disk213/xieqq/JINHUA138/Single_cell_analysis.16.CoVarNet")

####1. 在scRNA-seq数据中发现细胞模块####
#meta <- read.csv("meta.xls",row.names=1) 
meta <- read.csv("mymeta.csv",row.names=1) %>% column_to_rownames(var="cell_name")

rt_sp <- names(table(meta$sampleID))[table(meta$sampleID) >= 100]
meta <- meta[meta$sampleID %in% rt_sp, ]

mat_fq_raw <- freq_calculate(meta) 
mat_fq_norm <- freq_normalize(mat_fq_raw,normalize="minmax") #归一化

res <- nmf(mat_fq_norm, rank=2:10, method="nsNMF", seed=rep(123456, 6), .options="vp")
plot(res)
pdf("nmf_result.pdf", width=8, height=6) 
plot(res)
dev.off() 

rho <- setNames(res[["measures"]][["cophenetic"]], res[["measures"]][["rank"]])
find_pk <- function(rho) {
  ranks <- as.numeric(names(rho))
  for(i in 3:(length(rho)-1)) {
    if(rho[i-2] < rho[i-1] && rho[i-1] < rho[i] && rho[i] > rho[i+1]) {
      return(ranks[i])
    }
  }
  return(NA)
}
K <- find_pk(rho)
K=6
NMF_K <- nmf(mat_fq_norm, K, method="nsNMF", seed=rep(77, 6), nrun=30)
colnames(basis(NMF_K))=paste0("CM", sprintf("%02d", 1:K)) #basis(NMF_K)是占比
rownames(coef(NMF_K))=paste0("CM", sprintf("%02d", 1:K))
saveRDS(NMF_K,"NMF_K.rds")
#NMF_K <- readRDS("NMF_K.rds")

basis <- basis(NMF_K) %>% as.data.frame(.)
write.csv(basis,"basis.csv")
W_norm <- t(apply(basis, 1, function(x) x / sum(x)))
W_top10 <- apply(basis, 2, function(x) names(sort(x, decreasing=TRUE))[1:10])
write.csv(W_top10,"W_top10.csv")

coef <- coef(NMF_K) %>% as.data.frame(.)%>%t(.) #各模块在各样本中的活性

gr.weight_all(NMF_K) 
pdf("weight_all.pdf", width=4, height=6) 
gr.weight_all(NMF_K)
dev.off() 

gr.weight_top(NMF_K,num=15) 
pdf("weight_top15.pdf", width=10, height=6) 
gr.weight_top(NMF_K,num=15)
dev.off() 

h <- coef(NMF_K)
max_cm <- apply(h, 2, function(x) rownames(h)[which.max(x)])
max_cm <- gsub("CM", "CMT", max_cm)

gr.distribution(NMF_K,meta=meta,group="PCK1_exp")  
pdf("weight_distribution_subCluster.pdf", width=8, height=6) 
gr.distribution(NMF_K,meta=meta,group="subCluster") 
dev.off() 
pdf("weight_distribution_tissue.pdf", width=8, height=6) 
gr.distribution(NMF_K,meta=meta,group="tissue") 
dev.off() 

cor_pair <- pair_correlation(mat_fq_raw,method="pearson")
saveRDS(cor_pair,"cor_pair.rds")
filtered_cor_pair <- cor_pair %>% filter(correlation>=0.2,pval_fdr<0.05)
write.csv(filtered_cor_pair,"cor_pair_cor0.2_fdr0.05.csv")
#cor_pair <- readRDS("cor_pair.rds")


network <- cm_network(NMF_K,cor_pair,top_n=10,corr=0.2,fdr=0.05)
each <- network$each
gr.igraph_each(each,Layout=layout_in_circle) 
pdf("igraph_each.pdf", width=8, height=6) 
gr.igraph_each(each,Layout=layout_in_circle) 
dev.off() 

global <- network$global
par(plt=c(0,1,0,1),fig = c(0,1,0,1))
gr.igraph_global(global,Layout=layout_with_fr)
pdf("igraph_global.pdf", width=8, height=6) 
gr.igraph_global(global,Layout=layout_with_fr)
dev.off() 


head(network$filter)
saveRDS(network,"network.rds")
#network <- readRDS("network.rds")


####2. 分析####
library(dtw)
library(ggpubr)

coef <- scoef(NMF_K) %>% t(.) %>% as.data.frame(.) %>% rownames_to_column(var="sampleID")

PCK1_exp <- meta %>% select(sampleID,PCK1) %>%
  group_by(sampleID) %>%
  summarise(PCK1_mean=mean(PCK1, na.rm=TRUE)) %>%
  mutate(PCK1_group=case_when(PCK1_mean<0.035 ~ "Low", PCK1_mean>=0.035&PCK1_mean<0.218 ~ "Mid", PCK1_mean>=0.218 ~ "High")) %>% 
  mutate(PCK1_group=factor(PCK1_group, levels=c("Low","Mid","High")))

data <- data.frame(sampleID=meta$sampleID) %>% unique(.) %>% 
  left_join(.,PCK1_exp,by=c("sampleID")) %>% 
  left_join(.,coef,by=c("sampleID"))

cols <- c("Low"="#d1495b", "Mid"="#00798c", "High"="#edae49")
plot <- ggplot(data, aes(x=PCK1_group, y=CM03)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(aes(color=PCK1_group), width=0.15, size=2) +
  stat_compare_means(comparisons=list(c("Low","Mid"),c("Low","High"),c("Mid","High")),
                     method="t.test",label="p.format", hide.ns=TRUE, label.y=1) +
  scale_color_manual(values=cols) +
  labs(x="", y="CM03 activity") +
  theme_bw()+
  theme(axis.text.x=element_text(color="black",size=10),
        axis.text.y=element_text(color="black",size=10),
        axis.title.y=element_text(color="black",size=12),
        legend.text=element_text(color="black",size=10),
        legend.title=element_text(color="black",size=12),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
ggsave(filename="Boxplot_CM03_PCK1group.pdf", plot=plot, width=7, height=4)


data <- data.frame(sampleID=meta$sampleID,age=meta$age,tissue=meta$tissue) %>% unique(.) %>% 
  mutate(age=gsub("Adult","240",age)) %>%
  mutate(age=as.numeric(gsub("d","",age)),
         agegroup=case_when(age<30 ~ "G1", age>=30&age<60 ~ "G2", age>=60&age<90 ~ "G3", 
                            age>=90&age<180 ~ "G4", age>=180 ~ "G5")) %>%
  mutate(agegroup=factor(agegroup, levels=c("G1","G2","G3","G4","G5"))) %>% 
  left_join(.,coef,by=c("sampleID"))

cols <- c("#3B5F7F","#7AA6A1","#C4D7B2","#E8C39E","#D28B70")
plot <- ggplot(data, aes(x=agegroup, y=CM03)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(aes(color=agegroup), width=0.15, size=2) +
  stat_compare_means(comparisons=list(c("G1","G5"),c("G2","G5"),c("G3","G5"),c("G4","G5")),
                     method="t.test",label="p.format", hide.ns=TRUE, label.y=1) +
  scale_color_manual(values=cols) +
  labs(x="", y="CM03 activity") +
  theme_bw()+
  theme(axis.text.x=element_text(color="black",size=10),
        axis.text.y=element_text(color="black",size=10),
        axis.title.y=element_text(color="black",size=12),
        legend.text=element_text(color="black",size=10),
        legend.title=element_text(color="black",size=12),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
ggsave(filename="Boxplot_CM03_Agegroup.pdf", plot=plot, width=7, height=4)

