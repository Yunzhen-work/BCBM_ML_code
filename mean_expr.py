import pandas as pd
import numpy as np

expr_data = pd.read_csv(f"D:/Desktop/BCBM_by_ML/GSE125989.csv")
print(expr_data.head)

# 把Gene Symbol这列设为行名（索引）
expr_data_1 = expr_data.set_index('Gene Symbol')

# 对行名（索引）相同的行进行编组，并查看编组结果
grouped = expr_data_1.groupby(expr_data_1.index)
print(grouped.groups)

# 计算同组的平均值
avg_df = grouped.mean()

# 查看文件后发现第二列是无用数据，使用df.pop('列名')删掉
avg_df.pop('Unnamed: 0')
# avg_df.to_csv(f"D:/Desktop/BCBM_by_ML/test.csv")

###--------------------------------------------------------------------
# 处理带分隔符“ /// ”的基因，并复制整行到下一行，每个基因都复制一次
###--------------------------------------------------------------------
new_rows = []
num_gene = []
new_df = pd.DataFrame()
for row in avg_df.index:    
    # 如果行名里面带有分隔符' /// '，把基因按分隔符分割出来，并存入num_gene中
    if ' /// ' in row:
        num_gene = row.split(' /// ')
        # print(num_gene)
        
        # 逐步提取num_gene中的基因并复制后面的行，添加到新的dataframe中
        for i in num_gene:
            row_name = i
            row_data = avg_df.loc[row]
            # print(row_data.head)
            new_rows = pd.Series(row_data, name=row_name)
            new_df = pd.concat([new_df, pd.DataFrame([new_rows])]) 
            # print(new_df.head)

    # 如果这个行名里不包含分隔符（即只有一个基因），直接复制整行到新的dataframe中
    else:
        new_rows = pd.Series(avg_df.loc[row], name=row)
        new_df = pd.concat([new_df, pd.DataFrame([new_rows])])
   
# new_df.to_csv(f"D:/Desktop/BCBM_by_ML/test_1.csv")

###--------------------------------------------------------------------


# 处理完带分隔符的行后，表里又出现了重复的行，需要去重复
grouped_2 = new_df.groupby(new_df.index)
avg_df_2 = grouped_2.mean()
avg_df_2.to_csv(f"D:/Desktop/BCBM_by_ML/GSE125989_expr.csv")






