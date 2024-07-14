"""
在python中对两矩阵做相关性分析
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# 读入数据
exp_matrix_1 = pd.read_csv(f"D:/Desktop/BCBM_by_ML/test-1.csv")
exp_matrix_2 = pd.read_csv(f"D:/Desktop/BCBM_by_ML/test-2.csv")

# 设置行名(索引)
exp_1 = exp_matrix_1.set_index(exp_matrix_1.columns[0])
exp_2 = exp_matrix_2.set_index(exp_matrix_2.columns[0])
# print(exp_2.index)

# 开始进行相关性分析
new_rows = []
new_df = pd.DataFrame()
for i in exp_1.index:
    for j in exp_2.index:
        inputdata_1 = exp_1.loc[i]
        inputdata_2 = exp_2.loc[j]
        rho, p_value = spearmanr(inputdata_1, inputdata_2)
        # 把结果依次放入新的行里，一行行添加到新的Dataframe中
        new_data = {'gene_1':i, 'gene_2':j, 'SCC':rho, 'p_value':p_value}
        new_rows = pd.Series(new_data)
        new_df = pd.concat([new_df, pd.DataFrame([new_rows])]) 

# 由于前面没设置索引，这里第一列全显示0，我们把它删掉再输出
new_df_1 = new_df.reset_index(drop=True)
new_df_1.to_csv(f"D:/Desktop/BCBM_by_ML/test.csv")
# print(new_df)

# 如果需要对P值取小数点后几位数，可以用p_value:.4f
# print(f"斯皮尔曼相关系数: {rho:.4f}")
# print(f"p值: {p_value:.4f}")
### -------------------------------------------------------------------------
