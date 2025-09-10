
#### 细胞比例 ####

setwd("/disk213/xieqq/JINHUA138.sc/BayesPrism")

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)


fraction1 <- read.csv("./CellLineage_small/fraction_CellLineage_small.csv",row.names=1)
fraction2 <- read.csv("./CellLineage_large/fraction_CellLineage_large.csv",row.names=1)
fraction <- rbind(fraction1,fraction2)
fraction$B_T_lymphocytes <- fraction$Bcells+fraction$Plasma+fraction$Mesenchymal

fraction$Sample <- rownames(fraction)
fraction <- separate(fraction, col="Sample", into=c("breed","time","intestine"), sep="_", remove=T)

new_fraction <- fraction %>% 
  select("EECs","Tuft","Colonocytes","Goblet","BEST4enterocytes","TA","Stem","Progenitor","Enterocytes","time","intestine")
new_fraction <- fraction %>% 
  select("Epithelial","B_T_lymphocytes","T_ILC_NKcells","Myeloid","Endothelial","Neuronal","time","intestine")

data = new_fraction
row_sums <- rowSums(data[, 1:(ncol(data)-2)], na.rm = TRUE)
data[, 1:(ncol(data)-2)] <- data[, 1:(ncol(data)-2)] / row_sums

# 绘图
data <- melt(data, id.vars = c("time","intestine"), measure.vars = colnames(data)[1:(ncol(data)-2)])
data$time <- factor(data$time,levels=c("0","60","90","180","240"),labels=c("0d","60d","90d","180d","240d"))
data$intestine <- factor(data$intestine,levels=c("Small","Large"))

p <- ggplot(data, aes(x=time, y=value, fill=intestine)) +
  geom_boxplot(position=position_dodge(width=0.8), na.rm=TRUE) +
  scale_fill_manual(values=c("Large"="#B499C9", "Small"="#F0BDBC"))+
  facet_wrap(~ variable, scales="free_y") +
  labs(y="Proportion", x="time", fill="intestine") +
  theme_bw()+
  theme(axis.text.x=element_text(color="black", size=12),
        axis.text.y=element_text(color="black", size=12),
        axis.title.x=element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=14),
        legend.text=element_text(color="black",size=12),
        legend.title=element_text(color="black",size=12),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
p
ggsave(filename="fraction_CellLineage.pdf", plot=p, width=10, height=6)

summary_data <- fraction %>%
  group_by(time, intestine) %>%
  summarise(count = n()) %>%
  ungroup()
summary_data$time <- factor(summary_data$time,levels=c("0","60","90","180","240"),labels=c("0d","60d","90d","180d","240d"))
p2 <- ggplot(summary_data, aes(x = time, y = count, fill = intestine)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=c("Large"="#B499C9", "Small"="#F0BDBC"))+
  labs(title = "", x = "Time", y = "Count", fill = "Intestine Type") +
  theme_bw()+
  theme(axis.text.x=element_text(color="black", size=12),
        axis.text.y=element_text(color="black", size=12),
        axis.title.x=element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=14),
        legend.text=element_text(color="black",size=12),
        legend.title=element_text(color="black",size=12),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
p2
ggsave(filename="BayesPrism_sample_num.pdf", plot=p2, width=6, height=5)