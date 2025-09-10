
#### 2.tissue-specific genes ####
setwd("/disk213/xieqq/JINHUA138/Transcriptome_analysis.1.GTEx_Tissues")
sample <- list.files(path = "/disk191_1/zzy/GTExData/pigRNA",pattern="[.]ntpm.txt$", full.names=T)  #pareto scale
gene_list <-NULL
for (s in sample){
  Sample <- read.table(file=s,sep = ",",header=T)
  gene <- as.data.frame(Sample[,1])
  gene_list <- rbind(gene_list,gene)
  gene_list <- unique(gene_list)
}
colnames(gene_list) <- c('GENE_ID')
head(gene_list)

data <- gene_list
for (s in sample){
    a=gsub("/disk191_1/zzy/GTExData/pigRNA/","",s)
    tissue=gsub("[.]ntpm.txt","",a)
    Sample <- read.table(file=s,sep = ",",header=T)
    Sample$average <- rowMeans(Sample[,-1])
    data <- left_join(data,Sample[,c("GENE_ID","average")],by="GENE_ID")
    colnames(data)[ncol(data)] <- tissue
}
data <- arrange(data,GENE_ID)
write.csv(data,file = "GTExData_tissue_ntpm.csv",quote =F,row.names =F)

data=read.csv("GTExData_tissue_ntpm.csv")
cols_remain<- c("GENE_ID","Adipose","Muscle","Blood","Lymph_node","Large_intestine","Small_intestine",
                "Duodenum","Jejunum","Ileum","Colon",
                "Brain","Heart","Kidney","Liver","Lung","Spleen","Ovary","Uterus") 
data <- data[ ,colnames(data) %in% cols_remain]
write.csv(data,file = "GTExData_select_tissue_ntpm.csv",quote =F,row.names =F)

# Z
data$mean <- apply(data[,-1], 1, mean, na.rm = TRUE)
data$sd <- apply(data[,-1], 1, sd, na.rm = TRUE)
z_score <- (data[, -1]-data$mean)/data$sd #z-score
z_score <- cbind(data$GENE_ID,z_score)
colnames(z_score)[1] <- 'GENE_ID'
write.csv(z_score,file = "GTExData_select_tissue_ntpm_z-score.csv",quote =F,row.names =F)
z_score_abs <- abs((data[, -1]-data$mean)/data$sd) #|z-score|
z_score_abs <- cbind(data$GENE_ID,z_score_abs)
colnames(z_score_abs)[1] <- 'GENE_ID'
write.csv(z_score_abs,file = "GTExData_select_tissue_ntpm_|z-score|.csv",quote =F,row.names =F)

z_score = read.csv("GTExData_select_tissue_ntpm_z-score.csv")
z_score_abs = read.csv("GTExData_select_tissue_ntpm_|z-score|.csv")

max_vals <- apply(z_score_abs[,-1], 1, function(row) {
  max_val <- max(row, na.rm = TRUE)
  max_col <- names(row)[which.max(row)]
  c(max_val, max_col)
})
second_max_vals <- apply(z_score_abs[,-1], 1, function(row) {
  sorted_row <- sort(row, decreasing = TRUE)
  ifelse(length(sorted_row) > 1, sorted_row[2], NA)
})
tissue_specific_genes = as.data.frame(z_score_abs[,1])
tissue_specific_genes$max_val <- sapply(max_vals, "[", 1)
tissue_specific_genes$second_max_val <- second_max_vals
tissue_specific_genes$max_val_col <- sapply(max_vals, "[", 2)
tissue_specific_genes[, 4] <- ifelse(tissue_specific_genes[, 2] > 3*tissue_specific_genes[, 3], tissue_specific_genes[, 4], NA)
colnames(tissue_specific_genes) <- c('gene_id','nTPM_Zscore','second_max_val','tissue')

library(dplyr)
gene <- read.csv("/disk213/xieqq/JINHUA138/Genome_analysis.1.ROH/protein_coding_gene_JH.csv",header = T)
qtf.inp <- import_gff_gtf(db_file="/disk213/xieqq/Software/GALLO/Sus_scrofa.Sscrofa11.1.109.gtf",file_type="gtf")
tissue_specific_genes <- left_join(tissue_specific_genes,qtf.inp[,c(6,7)],by="gene_id")
tissue_specific_genes$gene_merge <- tissue_specific_genes$gene_name
tissue_specific_genes[is.na(tissue_specific_genes[,6]), 6] <- tissue_specific_genes[is.na(tissue_specific_genes[,6]), 1]
tissue_specific_genes=tissue_specific_genes[,c(1,6,4,2)]
tissue_specific_genes=arrange(tissue_specific_genes,tissue)

all_gene = tissue_specific_genes
write.csv(all_gene,file = "GTExData_select_tissue_ntpm_|z-score|_all_gene.csv",quote =F,row.names =F)
select_gene = tissue_specific_genes
select_gene <- subset(select_gene, gene_merge %in% gene$gene_merge)
write.csv(select_gene,file = "GTExData_select_tissue_ntpm_|z-score|_select_gene.csv",quote =F,row.names =F)

z_score_all_genes <- left_join(z_score,qtf.inp[,c(6,7)],by=c("GENE_ID"="gene_id"))
z_score_all_genes$gene_merge <- z_score_all_genes$gene_name
z_score_all_genes[is.na(z_score_all_genes[,ncol(z_score_all_genes)]), ncol(z_score_all_genes)] <- z_score_all_genes[is.na(z_score_all_genes[,ncol(z_score_all_genes)]), 1]
z_score_all_genes <- z_score_all_genes[,c(1,ncol(z_score_all_genes),2:(ncol(z_score_all_genes)-4))]
write.csv(z_score_all_genes,file = "GTExData_select_tissue_ntpm_z-score_all_gene_matrix.csv",quote =F,row.names =F)
matched_rows <- subset(z_score_all_genes, gene_merge %in% gene$gene_merge)
write.csv(matched_rows,file = "GTExData_select_tissue_ntpm_z-score_select_gene_matrix.csv",quote =F,row.names =F)

z_score_abs_all_genes <- left_join(z_score_abs,qtf.inp[,c(6,7)],by=c("GENE_ID"="gene_id"))
z_score_abs_all_genes$gene_merge <- z_score_abs_all_genes$gene_name
z_score_abs_all_genes[is.na(z_score_abs_all_genes[,ncol(z_score_abs_all_genes)]), ncol(z_score_abs_all_genes)] <- z_score_abs_all_genes[is.na(z_score_abs_all_genes[,ncol(z_score_abs_all_genes)]), 1]
z_score_abs_all_genes <- z_score_abs_all_genes[,c(1,ncol(z_score_abs_all_genes),2:(ncol(z_score_abs_all_genes)-4))]
write.csv(z_score_abs_all_genes,file = "GTExData_select_tissue_ntpm_|z-score|_all_gene_matrix.csv",quote =F,row.names =F)
matched_abs_rows <- subset(z_score_abs_all_genes, gene_merge %in% gene$gene_merge)
write.csv(matched_abs_rows,file = "GTExData_select_tissue_ntpm_|z-score|_select_gene_matrix.csv",quote =F,row.names =F)

#### nTPM pheatmap ####
setwd("/disk213/xieqq/JINHUA138/Transcriptome_analysis.1.GTEx_Tissues")
Sample <- read.csv("GTExData_select_tissue_ntpm.csv")
qtf.inp <- import_gff_gtf(db_file="/disk213/xieqq/Software/GALLO/Sus_scrofa.Sscrofa11.1.109.gtf",file_type="gtf")
Sample <- left_join(Sample,qtf.inp[,c(6,7)],by=c("GENE_ID"="gene_id"))
Sample$gene_merge <- Sample$gene_name
Sample[is.na(Sample[,ncol(Sample)]), ncol(Sample)] <- Sample[is.na(Sample[,ncol(Sample)]), 1]
write.csv(Sample,file = "GTExData_select_tissue_ntpm_all_gene_matrix.csv",quote =F,row.names =F)
gene <- read.csv("/disk213/xieqq/JINHUA138/Genome_analysis.1.ROH/protein_coding_gene_JH.csv",header = T)
Sample <- subset(Sample, gene_merge %in% gene$gene_merge)
write.csv(Sample,file = "GTExData_select_tissue_ntpm_select_gene_matrix.csv",quote =F,row.names =F)

# all gene
Sample <- read.csv("GTExData_select_tissue_ntpm_all_gene_matrix.csv")
df <- Sample[,-c(ncol(Sample),(ncol(Sample)-1))]
rownames(df)=df[,1]
df <- df[,-1]
df[is.na(df)]=0
df <- log10(df+1)

wss <- (nrow(df)-1)*sum(apply(df,2,var))
for (i in 2:500) 
  wss[i] <- sum(kmeans(df,centers=i)$withinss)
set.seed(999)
pdf(file = "kmeans_allgene.pdf")   # kmeans_k
plot(1:40, wss[1:40], type="b", xlab="Number of Clusters",ylab="Within groups sum of squares") 
dev.off()

pheatmap <- pheatmap(df, kmeans_k = 38, 
                     color = brewer.pal(8,"Blues")[1:8],
                     show_colnames = T,
                     #cutree_rows = 6,cutree_cols = 8,
                     angle_col=90)
pdf(file = "pheatmap_allgene.pdf",width=7, height=10)
pheatmap
dev.off()

kmeans=as.data.frame(rownames(df))
kmeans$cluster=NA
colnames(kmeans)[1]="gene"
for (i in rownames(df)){
  kmeans[which(kmeans$gene==i),2]=pheatmap[["kmeans"]][["cluster"]][[i]]
}
kmeans <- cbind(kmeans,Sample[,21])
colnames(kmeans)[3] <- "gene_merge"
write.csv(kmeans,file = "kmeans_allgene.csv",quote =F,row.names =F)


#### z score pheatmap ####
# all gene
Sample <- read.csv("GTExData_select_tissue_ntpm_z-score_all_gene_matrix.csv",row.names=c("GENE_ID"))
df <- Sample[,-1]
df[is.na(df)]=0

library(factoextra)
wss <- (nrow(df)-1)*sum(apply(df,2,var))
for (i in 2:100) 
  wss[i] <- sum(kmeans(df,centers=i)$withinss)
set.seed(999)
pdf(file = "kmeans_allgene.pdf")   
plot(1:40, wss[1:40], type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
dev.off()

library(RColorBrewer)
library(pheatmap)
pheatmap <- pheatmap(df, kmeans_k = 37,
                     color = brewer.pal(8,"Blues")[1:8],
                     show_colnames = T,
                     #cutree_rows = 6,cutree_cols = 8,
                     angle_col=90)
pdf(file = "pheatmap_allgene.pdf",width=7, height=10)
pheatmap
dev.off()

kmeans=as.data.frame(rownames(df))
kmeans$cluster=NA
colnames(kmeans)[1]="gene"
for (i in rownames(df)){
  kmeans[which(kmeans$gene==i),2]=pheatmap[["kmeans"]][["cluster"]][[i]]
}
kmeans <- cbind(kmeans,Sample[,1])
colnames(kmeans)[3] <- "gene_merge"
write.csv(kmeans,file = "kmeans_allgene.csv",quote =F,row.names =F)
