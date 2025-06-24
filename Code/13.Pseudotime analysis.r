install.packages("installr")
install.packages("stringr") 
library(stringr)
library(installr)
install.Rtools()

Sys.which("make")

install.packages('Cairo')

## make sure you remove the old version of TSCAN if you've downloaded it before
if ("TSCAN" %in% rownames(installed.packages())) remove.packages('TSCAN')
## download dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  BiocManager::install('ComplexHeatmap')
if (!require("devtools"))
  install.packages("devtools")
## downlaod the most updated TSCAN from Github
devtools::install_github("zji90/TSCAN") 
## downlaod the most updated Lamian from Github
devtools::install_github("Winnie09/Lamian")

install.packages ('Seurat')



#### Module 1: tree variability ####
options(warn=-1)
suppressMessages(library(Lamian)) # load in Lamian
library(ggplot2)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))

setwd("/disk213/xieqq/JINHUA138.sc/Lamian")

pbmc = readRDS("/disk213/xieqq/JINHUA138.sc/RDS/Epithelial.rds")
umap <- subset(pbmc, CellType %in% c("TA","Stem","Progenitor"))

plotdir <- "./tree_variability/"

cell_type_cols=c("#34a0a4","#184e77","#1a759f")
pdf(paste0(plotdir, 'pca.pdf'), width=7, height=5)
DimPlot(umap, reduction="umap", group.by="CellType", cols=cell_type_cols, label=TRUE, label.size=3, repel=TRUE)
dev.off()

pca <- as.matrix(umap@reductions$pca@cell.embeddings)
expression <- as.matrix(umap@assays$RNA@data)
cellanno <- data.frame(cell=rownames(umap@meta.data), sample=umap@meta.data$SEGMENT.TIME, celltype=umap@meta.data$CellType)
rownames(cellanno)=cellanno$cell

### infer tree structure
res = infer_tree_structure(pca = pca,
                           expression = expression,
                           cellanno = cellanno,
                           origin.marker = c('LGR5'),
                           origin.celltype = 'Stem',
                           number.cluster = 6,
                           plotdir = plotdir,
                           xlab='Principal component 1',
                           ylab='Principal component 2')

pdf(paste0(plotdir, 'tree_structure.pdf'), width=6,height=5)
plotmclust(res, cell_point_size=0.5, x.lab='Principal component 1', y.lab = 'Principal component 2')
dev.off()

pseudotime0 <- as.data.frame(as.matrix(res[["pseudotime"]]))
pseudotime0$V2 <- names(res[["pseudotime"]])
pseudotime0 <- pseudotime0[!duplicated(pseudotime0$V2),]
pseudotime <- pseudotime0$V1
names(pseudotime) <- pseudotime0$V2
write.csv(pseudotime,paste0(plotdir, "pseudotime.csv"))

### evaluate branch uncertainty
result <- evaluate_uncertainty(res, n.permute=10)
saveRDS(result, paste0(plotdir, 'result.rds'))


#### Module 3: Trajectory differential tests about gene expression ####
options(warn=-1)
suppressMessages(library(Lamian)) # load in Lamian
library(ggplot2)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
library(cluster)
library(factoextra)
library(dplyr)
library(pheatmap)
library(gridExtra)


setwd("/disk213/xieqq/JINHUA138.sc/Lamian")

pbmc = readRDS("/disk213/xieqq/JINHUA138.sc/RDS/Epithelial.rds")
pbmc <- subset(pbmc, CellType %in% c("Colonocytes","Enterocytes"))
unique(pbmc$SEGMENT)

plotdir <- "./trajectory_tests1/"
segment <- "cecum"
umap <- subset(pbmc, SEGMENT %in% segment)

pdf(paste0(plotdir, segment, '_pca.pdf'), width=6, height=5)
DimPlot(umap, reduction="pca", group.by="SEGMENT.TIME", label=TRUE, label.size=3, repel=TRUE)
dev.off()

pca <- as.matrix(umap@reductions$pca@cell.embeddings)
expression <- as.matrix(umap@assays$RNA@data)
cellanno <- data.frame(cell=rownames(umap@meta.data), sample=umap@meta.data$SEGMENT.TIME, 
                       celltype=sapply(strsplit(umap@meta.data$SEGMENT.TIME, "-"), function(x) x[2]))
rownames(cellanno)=cellanno$cell

res = infer_tree_structure(pca = pca, expression = expression, cellanno = cellanno, 
                           origin.marker=c('ANPEP','FABP2'), origin.celltype=c("0"),
                           number.cluster = 5, plotdir = paste0(plotdir,segment,"_"),
                           xlab='Principal component 1', ylab='Principal component 2')

pdf(paste0(plotdir,segment,"_tree_structure.pdf"), width=5,height=5)
plotmclust(res,x=1,y=2,cell_point_size=0.5)
dev.off()

pseudotime <- res[["pseudotime"]]
write.csv(pseudotime,paste0(plotdir,segment, "_pseudotime.csv"))

design = data.frame(intercept=1,group=unique(cellanno$sample))
rownames(design) <- design$group
design$group <- sapply(strsplit(design$group, "-"), function(x) x[2])

saveh5(expr=expression, pseudotime=pseudotime, cellanno=cellanno, path=paste0(plotdir,segment,"_trajectory.h5"))

Res <- lamian_test(expr=expression, cellanno=cellanno, pseudotime=pseudotime, design=design, 
                   test.type='Time', test.method="permutation", testvar=2, permuiter=5, ncores=1)

## determine the TDE genes as the genes with fdr.overall <0.05
diff_gene <- Res$statistics[Res$statistics[, 1] < 0.05,]
diffgene <- diff_gene %>% rownames()

## population fit
num.timepoint=max(pseudotime)
Res$populationFit <- getPopulationFit(Res, gene=diffgene, type='time', num.timepoint=num.timepoint)
## clustering
Res$cluster <- clusterGene(Res, gene=diffgene, type='time', k.auto=T, method="kmeans")
max(Res$cluster)
pdf(paste0(plotdir,segment,'_cluster_mean.pdf'), width=5, height=6)
plotClusterMean(Res, cluster=Res$cluster, type='time')
dev.off()

fit <- testTDEHm_i(Res, showRowName=F, subsampleCell=F, showCluster=T, type='time')
colnames(fit) <- c(1:num.timepoint)
write.csv(fit, paste0(plotdir, segment,"_fit.pseudotime.csv"))

outgene <- data.frame(gene=rownames(fit),order=c(1:nrow(fit)))
outcluster <- data.frame(cluster=Res$cluster)
outcluster$gene <- rownames(outcluster)
diff_gene$gene <- rownames(diff_gene)
outgene <- left_join(outgene,outcluster,by="gene")
outgene <- left_join(outgene,diff_gene,by="gene")
write.csv(outgene, paste0(plotdir, segment,"_gene_fdr.csv"))

## save png
col.expression = brewer.pal(n = 8, name = "Pastel1")[seq_len(2)]
names(col.expression) = c('Original', 'Model Fitted')
colann.fit <- data.frame(pseudotime = colnames(fit), expression = 'Model Fitted', stringsAsFactors = FALSE)
col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(num.timepoint)
names(col.pseudotime) = colnames(fit)
clu <- Res$cluster
if (length(unique(clu)) < 8) {
  col.clu = brewer.pal(8, 'Set1')[seq_len(length(unique(clu)))]
} else {
  col.clu = colorRampPalette(brewer.pal(8, 'Set1'))[seq_len(length(unique(clu)))]
}
names(col.clu) = unique(clu)
annotation_colors = list(pseudotime = col.pseudotime,expression = col.expression,cluster = col.clu)
cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
rowann = data.frame(cluster = as.character(clu),stringsAsFactors = FALSE)
rownames(rowann) = names(clu)
rowann <- rowann[rownames(fit), , drop = FALSE]
cellWidthTotal = 250;cellHeightTotal = 400

pheatmap(fit, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         color = cpl, annotation_col = colann.fit, annotation_row = rowann, annotation_colors = annotation_colors,
         cellwidth = cellWidthTotal / ncol(fit), cellheight = cellHeightTotal / nrow(fit), border_color = NA,
         silent = TRUE, filename = paste0(plotdir,segment,'_pheatmap.png'))

uniform_samples_seq <- ceiling(seq(from = 1, to = ncol(fit), length.out = 1000))
pheatmap(fit[,uniform_samples_seq], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         color = cpl, annotation_row = rowann, annotation_colors = annotation_colors,
         cellwidth = cellWidthTotal / length(uniform_samples_seq), cellheight = cellHeightTotal / nrow(fit), border_color = NA,
         silent = TRUE, filename = paste0(plotdir,segment,'_pheatmap_f.png'))

#### cluster ####
file_list <- list.files(path="/disk213/xieqq/JINHUA138.sc/Lamian/trajectory_tests", pattern="fdr.csv", full.names=TRUE)
pbmc <- readRDS("/disk213/xieqq/JINHUA138.sc/RDS/Epithelial.rds")
for (file in file_list) {
  segment <- gsub(".*/|(\\_gene_fdr.csv$)", "", file)
  seurat_obj <- subset(pbmc, SEGMENT %in% segment)
  data <- read.csv(file,check.names=F,row.names=1)
  FetchData <- data.frame(cells=rownames(seurat_obj@meta.data))
  for (i in unique(data$cluster)){
    genes_list <- data$gene[which(data$cluster==i)]
    means <- data.frame(cells=colnames(seurat_obj[["RNA"]]@data[genes_list, ]),means=colMeans(seurat_obj[["RNA"]]@data[genes_list, ]))
    rownames(means) <- NULL
    FetchData <- left_join(FetchData,means,by="cells")
  }
  colnames(FetchData) <- c("cells",unique(data$cluster))
  MEs <- FetchData[,-1]
  rownames(MEs) <- rownames(FetchData)
  mods <- colnames(MEs)
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  FetchData$SEGMENT.TIME = seurat_obj@meta.data$SEGMENT.TIME
  FetchData$CellType = seurat_obj@meta.data$CellType
  write.csv(FetchData,paste0("/disk213/xieqq/JINHUA138.sc/Lamian/trajectory_tests/FetchData/Cluster_FetchData_",segment,".csv"))
  seurat_obj@meta.data$TIME <- factor(seurat_obj@meta.data$TIME, levels=c("240","180","90","60","0")) #倒序
  p1 <- DotPlot(seurat_obj, features=mods, group.by="TIME")+RotatedAxis()+scale_color_gradient2(high="#08519C", low="#C6DBEF")
  ggsave(filename=paste0("/disk213/xieqq/JINHUA138.sc/Lamian/trajectory_tests/FetchData/Cluster_TIME_",segment,".pdf"), plot=p1, width=8, height=8)
}
