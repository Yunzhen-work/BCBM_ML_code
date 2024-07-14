"""
Immune infiltration analysis for different groups by Cibersort
"""
# remove(list = ls()) # clear

# 加载包
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(preprocessCore)
library(furrr)
# library(Rfast)
library(e1071)

# 设置工作目录
setwd("D:/Desktop/BCBM_by_ML")

# 作者更新了CIBERSORT，这里用的是github上下载的2020年最新版
source("CIBERSORT.R")
source("CoreAlg.R")
source("doPerm.R")
source("untils.R")
source("zzz.R")

# 基因名那列设为 "Gene symbol"
# 行名不能有重复的，先处理一下
table1 = read.table("GSE125989_expr.txt", sep="\t", header = FALSE)
class(table1)
table1_unique <- table1[!duplicated(table1[, 1]), ]
write.table(table1_unique, file = "GSE125989_expr_1.txt", sep = "\t", quote = F, row.names = F, col.names = F) 

# 加载数据
exp_file <- "GSE125989_expr_1.txt"
LM22_file <- "LM22.txt"

# perm是置换次数，一般文章会选1000次，QN如果是芯片数据T，测序数据F
TME_temp = cibersort(LM22_file, exp_file, perm = 1000, QN = TRUE) 
TME_results = data.frame(GSE_sample = rownames(TME_temp), TME_temp) #列名错位，修正
write.table(TME_results, file = "TME_results.txt", sep = "\t", quote = F, row.names = F, col.names = T) 
# 给行名添加一个新的名称，这样输出文件就是对齐的


### ------------------------------------------------------------------------

### ------------------------------------------------------------------------
# 免疫浸润结果可视化
### ------------------------------------------------------------------------

# 做分组信息
phenotype = read.csv("sample_metadata.csv", sep=",", header=T, row.names=1)
group_list = phenotype$condition
table(group_list)

# 绘图数据预处理
TME_data = as.data.frame(TME_results[, 1:22]) #后三列不要

TME_data$group = group_list
TME_data$sample = row.names(TME_data)

TME_New = melt(TME_data) #融合数据
# head(TME_New)
TME_New = TME_New[,-1]
colnames(TME_New) = c("Group", "Sample", "Celltype", "Composition")

# 按免疫细胞占比中位数排序绘图
plot_order = TME_NeW[TME_New$Group=="Tumor",]%>%group_by(celltype)%>%
summarise(m=median(Composition))%>%arrange(desc(m))%>%pu1l(celltype)

TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)

# 绘制箱线图
if(T){
  mytheme = theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),axis,title =element_text(size = 12,color ="black"),
                 axis.text=element_text(size= 12,color ="black"),
                 panel.grid.minor.y=element_blank(),
                 panel.grid.minor.x=element_blank(),
                 axis.text.x=element_text(angle=45,hjust =1),panel.grid=element_blank(),
                 legend.position ="top",
legend.text=element_text(size= 12),
legend.title=element_text(size= 12)
)}

boX_TME = ggplot(TME_New,aes(x=Celltype,y= Composition))+
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+
  scale_fill_manual(values=c("#1CB4B8", "#EB7369"))+
  theme_classic()+mytheme +
  stat_compare_means(aes(group=Group),
                     label ="p.signif",
                     method ="wilcox.test",
                     hide.ns = T)
  
boX_TME
ggsave("D:/Desktop/BCBM_by_ML/BCBM_TME.pdf", boX_TME, height = 15, width = 25, units = "cm")                                                                                                                                        

### -----------------------------------------------------------------------------



