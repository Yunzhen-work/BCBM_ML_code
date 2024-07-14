"""
在R中，两个处理好的表达谱矩阵，两个矩阵之间的基因两两做相关性计算
"""

table1=read.table("gene_matrix-1.txt",header=TRUE)
table2=read.table("gene_matrix-2.txt",header=TRUE)

table1_NameColumn = as.vector(table1[, 1])
# 用-1取矩阵中的全部表达值
table1=table1[,-1]

table2_NameColumn = as.vector(table2[, 1])
table2=table2[,-1]

n1=nrow(table1)
n2=nrow(table2)

# 新建一个matrix装结果
output <- matrix(0,n1*n2,5)
colnames(output) <- c("gene1","gene2","pearson_value","p_value","FDR")

row=0
count=0
for(i in 1:n1){
for(j in 1:n2){
	x <- vector()
	x <- as.numeric(table1[i,])
	
	y <- vector()
	y <- as.numeric(table2[j,])
	
	# 开始计算相关性
	cor <- cor.test(x, y, method = "pearson")
	p_value <- cor$p.value
	pearson_value <- as.numeric(cor$estimate)
	count=count+1
	output[count, 1] <- table1_NameColumn[i]
	output[count, 2] <- table2_NameColumn[j]
	output[count, 3] <- pearson_value
	output[count, 4] <- p_value
}	
}

output[,5] <- p.adjust(output[,4],method="BH")
write.table(output, "PCC.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
