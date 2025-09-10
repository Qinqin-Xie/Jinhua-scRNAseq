#### scProportionTest ####

library(scProportionTest)
library(Seurat)
library(ggplot2)

seurat_data <- readRDS("/disk213/xieqq/JINHUA138.sc/RDS/CellLineage.rds")

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
unique(seurat_data$IT)
order <- c("small_0","small_60","small_90","small_180","small_240",
           "large_0","large_60","large_90","large_180","large_240")
for (i in 1:(length(order)-1)){
  for (j in (i+1):length(order)){
    test_result <- permutation_test(prop_test,cluster_identity="CellLineage",sample_identity="IT",
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
data$result <- "Not_Significant"
data$result[which(data$FDR<0.05&abs(data$obs_log2FD)>1.5)] <- "Significant"
write.csv(data,"permutation_test_CellLineage.csv",row.names=F)
