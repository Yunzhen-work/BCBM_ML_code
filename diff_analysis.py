"""
差异表达分析及可视化
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

# 先把数据导入并处理好格式
avg_df_2 = pd.read_csv(f"D:/Desktop/BCBM_by_ML/GSE125989_expr.csv")
# print(avg_df_2.iloc[:,0]) #查看第一列
# 把第一列（基因名）设为index，由于第一列没有列名，用df.set_index(df.columns[0])
gene_exp = avg_df_2.set_index(avg_df_2.columns[0])
sample_metadata = pd.read_csv(f"D:/Desktop/BCBM_by_ML/sample_metadata.csv")
group_inf = sample_metadata.set_index(sample_metadata.columns[0])
# print(group_inf)

# 用样本分组信息来做分组矩阵
a = group_inf['condition'].tolist()
list1 = []
list2 = []
for i in a:
    if i == "brain metastases":
        x = 1
        y = 0
    else:
        x = 0
        y = 1
    list1.append(x)
    list2.append(y)
design = pd.DataFrame({'metastases': list1, 'primary': list2}, index=a)

# 构建差异比较矩阵
limma = importr("limma")
contrast_matrix = limma.makeContrasts("metastases-primary", levels = design)

# 调用R中的limma包进行差异表达分析
fit = limma.lmFit(gene_exp, design)
fit2 = limma.contrasts_fit(fit, contrast_matrix)
fit2 = limma.eBayes(fit2)
DEG_ot = limma.topTable(fit2, adjust="fdr", coef=1, n=float("inf"))
# print(DEG_ot)

# 因为前面用的limma包是R的，因此这里我们需要将结果转换回python对象
# 使用ro.conversion将结果转换回python对象
diff_genes = ro.conversion.rpy2py(DEG_ot)
# print(diff_genes)
# diff_genes.to_csv(f"D:/Desktop/BCBM_by_ML/GSE125989_diff_all.csv")

# 进一步提取差异表达基因
newlist1 = []
newlist2 = []
filtered_gene = pd.DataFrame()
Not_significant = pd.DataFrame()
# 把显著的结果放在diff_genes_filtered里，不显著的放在Not_significant里
for j in range(0, len(diff_genes.index)):
        if (abs(diff_genes['logFC'][j])>=1) and (diff_genes['adj.P.Val'][j]<=0.05):
             name1 = diff_genes.index[j]
             data_1 = diff_genes.iloc[j]
             newlist1 = pd.Series(data_1, name=name1)
             # print(newlist1)
             filtered_gene = pd.concat([filtered_gene, pd.DataFrame([newlist1])]) 
        else:
             name2 = diff_genes.index[j]
             data_2= diff_genes.iloc[j]
             newlist2 = pd.Series(data_2, name=name2)
             Not_significant = pd.concat([Not_significant, pd.DataFrame([newlist2])])

# 查看显著和不显著分别有多少
print(len(filtered_gene.index))
print(len(Not_significant.index))             

# 把筛选后的差异基因输出来             
# filtered_gene.to_csv(f"D:/Desktop/BCBM_by_ML/GSE125989_diff_filtered.csv")

# 上下调基因
newlist3 = []
up_list = pd.DataFrame()
for i in range(0, len(diff_genes.index)):
        if (diff_genes['logFC'][i]>=1) and (diff_genes['adj.P.Val'][i]<=0.05):
             name3 = diff_genes.index[i]
             data_3 = diff_genes.iloc[i]
             newlist3 = pd.Series(data_3, name=name3)
             up_list = pd.concat([up_list, pd.DataFrame([newlist3])])

newlist4 = []
down_list = pd.DataFrame()
for k in range(0, len(diff_genes.index)):
        if (diff_genes['logFC'][k]<-1) and (diff_genes['adj.P.Val'][k]<=0.05):
             name4 = diff_genes.index[k]
             data_4 = diff_genes.iloc[k]
             newlist4 = pd.Series(data_4, name=name4)
             down_list = pd.concat([down_list, pd.DataFrame([newlist4])])
             print(down_list)

print("上调基因数量：", len(up_list))
print("下调基因数量：", len(down_list))
print("差异表达基因数量：", len(up_list)+len(down_list))
             
# 把上下调基因输出来             
# up_list.to_csv(f"D:/Desktop/BCBM_by_ML/GSE125989_diff_up.csv")
# down_list.to_csv(f"D:/Desktop/BCBM_by_ML/GSE125989_diff_down.csv")

# 差异表达分析可视化——vocalno
import seaborn as sns
diff_genes['adj.P.Val_log'] = -np.log10(diff_genes['adj.P.Val'])
diff_genes['sig'] = 'normal'
diff_genes.loc[(diff_genes['logFC']>=1)&(diff_genes['adj.P.Val']<0.05),'sig'] = 'up'
diff_genes.loc[(diff_genes['logFC']<-1)&(diff_genes['adj.P.Val']<0.05),'sig'] = 'down'

plt.figure(figsize=(12, 8))
colors = ["#01c5c4","#ff414d","#686d76"]
sns.set_palette(sns.color_palette(colors))

sns.scatterplot(x='logFC', y='adj.P.Val_log', data=diff_genes,
                     hue='sig',
                     edgecolor=None, #node edgecolor
                     s=8, # size of node
                     )
plt.title('Vocalno Plot')
plt.xlabel('log2 Fold Change')
plt.ylabel('adj.P.Val_log')
# plt.show()







