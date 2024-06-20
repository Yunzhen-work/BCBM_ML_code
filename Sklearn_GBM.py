"""
使用sklearn，Gradient Boosting梯度提升(GBM)，留一法对样本进行二分类
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler # 做标准化
from sklearn import preprocessing
# from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_curve, auc, recall_score
import matplotlib.pyplot as plt


# 读入数据
brain_8gene = pd.read_csv(f"D:/Desktop/BCBM_by_ML/file_input/brain_73_8gene_38sample.csv")
metadata = pd.read_csv(f"D:/Desktop/BCBM_by_ML/file_input/sample_metadata_73_38sample.csv")

# ### -------------------------------------------------------------------------
# 把数据格式修理一下
# ### -------------------------------------------------------------------------
brain_8gene_T = brain_8gene.T #转置
new_columns = brain_8gene_T.iloc[0] # 第一行做列名
brain_8gene_T = brain_8gene_T.rename(columns=new_columns)
# # 删掉原来的第一行
brain_8gene_T = brain_8gene_T.drop(brain_8gene_T.index[0]) 
# print(brain_8gene_T)

metadata = metadata.set_index('Sample')
# ### -------------------------------------------------------------------------

### ---------------------------------------------------------------------------
# 开始做GBM
### ---------------------------------------------------------------------------

# 将特征df和标签df组合成X和y
X = brain_8gene_T.values
y = metadata.values

# 需要把标签改成（0，1），或者（-1，1）
# 这里需要用LabelEncoder()函数将响应变量的连续值简单地转换为分类值
lab = preprocessing.LabelEncoder()
y = lab.fit_transform(y)
y = 2*y -1
print(y)

# GBM+LOOVC
loo = LeaveOneOut()
loo.get_n_splits(X)
predictlable_list = []
reallable_list = []
scaler = StandardScaler()
X = scaler.fit_transform(X)
clf = GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1, random_state=0)
count_right_lable = 0
count = 0 #循环次数

Y_score_list = []
true_lables = []
pred_probs = []
for tran_index, test_index in loo.split(X):
    X_train, X_test = X[tran_index], X[test_index]
    Y_train, Y_test = y[tran_index], y[test_index]
    clf.fit(X_train, Y_train)
    predictlable_list.append(list(clf.predict(X_test)))
    reallable_list.append(list(Y_test))

    # 用来画后面的ROC
    Y_pred_prob = clf.predict_proba(X_test)[:,1]
    true_lables.extend(Y_test)
    pred_probs.extend(Y_pred_prob)

    if Y_test ==clf.predict(X_test):
        count_right_lable += 1
    count += 1
    print('第{}次循环'.format(count))
accuracy = count_right_lable/len(y)
print('-----循环结束------')
print('准确率为：%.2f%%'%(accuracy*100))
print('-----运行结束------')

### -----------------------------------------------------------------------------
# 计算F1_score, Accuracy, Sensitivity, Specificity等
### -----------------------------------------------------------------------------

# 辅助函数1：将一个任意嵌套的列表整理为一个列表
def flatten(nested):
    try:
        for sublist in nested:
            for element in flatten(sublist):
                yield element
    except TypeError:
        yield nested

# 辅助函数2：获取预测结果的TP TN FP FN 
def get_TpTN_FpFn(list1, list2):
    # list1 为真实的lable, list2为预测的lable
    reallable_list = list(flatten(list1))
    predictlable_list = list(flatten(list2))
    TP_count = 0
    TN_count = 0
    FP_count = 0
    FN_count = 0
    for i in range(len(reallable_list)):
        if reallable_list[i] == -1 and predictlable_list[i] == -1:
            TP_count += 1
        if reallable_list[i] == 1 and predictlable_list[i] == 1:
            TN_count += 1
        if reallable_list[i] == 1 and predictlable_list[i] == -1:
            FP_count += 1
        if reallable_list[i] == -1 and predictlable_list[i] == 1:
            FN_count += 1
    return TP_count, TN_count, FP_count, FN_count

# 开始计算    
TP_count, TN_count, FP_count, FN_count = get_TpTN_FpFn(reallable_list, predictlable_list)
print(TP_count, TN_count, FP_count, FN_count)
F1_score = (2 * TP_count) / (2 * TP_count + FP_count + FN_count)
ACC = (TP_count + FP_count) / (TP_count + FN_count + TN_count + FP_count)
SEN = TP_count / (TP_count + FN_count)
SPE = TN_count / (TN_count + FP_count)
PRE = TP_count / (TP_count + FP_count)

print('F1_score为：%.2f%%' % (F1_score * 100))
print('Accuracy为：%.2f%%' % (ACC * 100))
print('Sensitivity为：%.2f%%' % (SEN * 100))
print('Specificity为：%.2f%%' % (SPE * 100))
print('Precision为：%.2f%%' % (PRE * 100))
### ---------------------------------------------------------------------------

### ---------------------------------------------------------------------------
#  计算AUC，绘制ROC
### ---------------------------------------------------------------------------

fpr, tpr, threshold = roc_curve(true_lables, pred_probs)
roc_auc = auc(fpr, tpr)
plt.figure()
lw = 3

plt.plot(fpr, tpr, color='darkorange', lw=lw, label='Roc curve(area=%0.2f)'% roc_auc)
plt.plot([0,1], [0,1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()
### -----------------------------------------------------------------------------
