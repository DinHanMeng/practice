


#引用包
library(limma)
library(ggpubr)

expFile="GES77298.txt"      #表达数据文件
geneFile="interGenes.txt"     #基因列表文件
conFile="s1.txt"              #对照组的样品信息文件
treatFile="s2.txt"            #实验组的样品信息文件
setwd("")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
rt=avereps(data)

#如果数据没有取log2,会对数据自动取log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	rt[rt<0]=0
	rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

#读取样品信息文件，判断样品的类型
con=read.table(conFile, header=F, sep="\t", check.names=F)
treat=read.table(treatFile, header=F, sep="\t", check.names=F)
conData=data[,as.vector(con[,1])]
treatData=data[,as.vector(treat[,1])]
data=cbind(conData, treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#读取基因列表文件,提取疾病特征基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=F]

#设置比较组
Type=c(rep("Control",conNum), rep("Treat",treatNum))
my_comparisons=list()
my_comparisons[[1]]=levels(factor(Type))

#对基因进行循环，绘制差异箱线图
for(i in row.names(data)){
	#data[i,][data[i,]>quantile(data[i,], 0.99)]=quantile(data[i,], 0.99)
	rt1=data.frame(expression=data[i,], Type=Type)

	#对差异基因进行可视化，绘制箱线图
	boxplot=ggboxplot(rt1, x="Type", y="expression", color="Type",
				      xlab="",
				      ylab=paste(i, "expression"),
				      legend.title="",
				      palette = c("blue", "red"),
				      add = "jitter")+ 
		#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
		stat_compare_means(comparisons = my_comparisons)
		
	#输出图片
	pdf(file=paste0("boxplot.",i,".pdf"), width=5, height=4.5)
	print(boxplot)
	dev.off()
}
		



