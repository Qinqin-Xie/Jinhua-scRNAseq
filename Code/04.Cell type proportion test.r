setwd("/disk213/xieqq/JINHUA138/Single_cell_analysis.1.Cell_fraction")

#### scProportionTest ####

library(scProportionTest)
library(Seurat)
library(ggplot2)

#seurat_data <- system.file("extdata", "example_data.RDS", package = "scProportionTest")
#seurat_data <- readRDS(seurat_data)

#seurat_data <- readRDS("/disk213/xieqq/JINHUA138.sc/RDS/CellLineage.rds")
seurat_data <- readRDS("/disk213/xieqq/JINHUA138.sc/RDS/Epithelial.rds")

prop_test <- sc_utils(seurat_data)

test_result <- permutation_test(prop_test,
                                cluster_identity="CellType",
                                sample_1="du60", 
                                sample_2="du0",
                                sample_identity="BATCH") 


p <- permutation_plot(test_result, FDR_threshold=0.05, log2FD_threshold=2.5, order_clusters=F)
p[["data"]][["obs_log2FD"]]
p[["data"]][["FDR"]]
p[["data"]][["clusters"]]

#start
data <- as.data.frame(matrix(data=NA,nrow=0,ncol=7,dimnames=list(NULL,c("clusters","obs_log2FD","boot_CI_low","boot_CI_high","FDR","group1","group2"))))
unique(seurat_data$BATCH)
order <- c("du0","du60","du90","du180","du240","je0","je60","je90","je180","je240","il0","il60","il90","il180","il240",
           "ce0","ce60","ce90","ce180","ce240","co0","co60","co90","co180","co240")
for (i in 1:(length(order)-1)){
  for (j in (i+1):length(order)){
    test_result <- permutation_test(prop_test,cluster_identity="CellType",sample_identity="BATCH",
                                    sample_1=order[i],sample_2=order[j])
    p <- permutation_plot(test_result, order_clusters=F)
    newdata <- data.frame(clusters=p[["data"]][["clusters"]],
                          obs_log2FD=p[["data"]][["obs_log2FD"]],
                          boot_CI_low=p[["data"]][["boot_CI_2.5"]],
                          boot_CI_high=p[["data"]][["boot_CI_97.5"]],
                          FDR=p[["data"]][["FDR"]])
    newdata$group1 <- order[i]
    newdata$group2 <- order[j]
    data <- rbind(data,newdata)
  }
}
data$compare <- "BATCH"
data$result <- "Not_Significant"
data$result[which(data$FDR<0.05&abs(data$obs_log2FD)>1.5)] <- "Significant"
write.csv(data,"permutation_test_Epithelial_Batch.csv",row.names=F)

data <- as.data.frame(matrix(data=NA,nrow=0,ncol=7,dimnames=list(NULL,c("clusters","obs_log2FD","boot_CI_low","boot_CI_high","FDR","group1","group2"))))
unique(seurat_data$TIME)
order <- c("0","60","90","180","240")
for (i in 1:(length(order)-1)){
  for (j in (i+1):length(order)){
    test_result <- permutation_test(prop_test,cluster_identity="CellType",sample_identity="TIME",
                                    sample_1=order[i],sample_2=order[j])
    p <- permutation_plot(test_result, order_clusters=F)
    newdata <- data.frame(clusters=p[["data"]][["clusters"]],
                          obs_log2FD=p[["data"]][["obs_log2FD"]],
                          boot_CI_low=p[["data"]][["boot_CI_2.5"]],
                          boot_CI_high=p[["data"]][["boot_CI_97.5"]],
                          FDR=p[["data"]][["FDR"]])
    newdata$group1 <- order[i]
    newdata$group2 <- order[j]
    data <- rbind(data,newdata)
  }
}
data$compare <- "TIME"
data$result <- "Not_Significant"
data$result[which(data$FDR<0.05&abs(data$obs_log2FD)>2.5)] <- "Significant"
write.csv(data,"permutation_test_Epithelial_Time.csv",row.names=F)

data <- as.data.frame(matrix(data=NA,nrow=0,ncol=7,dimnames=list(NULL,c("clusters","obs_log2FD","boot_CI_low","boot_CI_high","FDR","group1","group2"))))
unique(seurat_data$SEGMENT)
order <- c("duodenum","jejunum","ileum","cecum","colon")
for (i in 1:(length(order)-1)){
  for (j in (i+1):length(order)){
    test_result <- permutation_test(prop_test,cluster_identity="CellType",sample_identity="SEGMENT",
                                    sample_1=order[i],sample_2=order[j])
    p <- permutation_plot(test_result, order_clusters=F)
    newdata <- data.frame(clusters=p[["data"]][["clusters"]],
                          obs_log2FD=p[["data"]][["obs_log2FD"]],
                          boot_CI_low=p[["data"]][["boot_CI_2.5"]],
                          boot_CI_high=p[["data"]][["boot_CI_97.5"]],
                          FDR=p[["data"]][["FDR"]])
    newdata$group1 <- order[i]
    newdata$group2 <- order[j]
    data <- rbind(data,newdata)
  }
}
data$compare <- "SEGMENT"
data$result <- "Not_Significant"
data$result[which(data$FDR<0.05&abs(data$obs_log2FD)>2.5)] <- "Significant"
write.csv(data,"permutation_test_Epithelial_Segment.csv",row.names=F)