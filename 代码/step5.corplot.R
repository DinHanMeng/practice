

#引用包
options(stringsAsFactors=F)
library(corrplot)
library(circlize)
library(limma)
library(PerformanceAnalytics)
#输入文件
setwd("")
expFile="diffGeneExp.txt"     #表达矩阵，去除正常控制
hub="gene.txt"           #核心基因
        
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

#提取目的基因
VENN=read.table(hub,sep = "\t",header = F,check.names = F)[,1]
rt1=rt[c(VENN),]
rt1=cbind(ID=rownames(rt1),rt1)
#保存目的基因表达
write.table(x = rt1,file = "AFexp.txt",quote = F,sep = "\t",col.names = T,row.names = F)

#读取输入文件
rt=read.table("AFexp.txt",sep="\t",header=T,check.names=F,row.names=1)
rt=t(rt)  #转置
M=cor(rt)     #相关型矩阵
res <- cor.mtest(rt)

#绘制相关性图形
pdf(file="corpot1.pdf",width=7,height=7)
corrplot(M,
         method = "circle",
         order = "hclust", #聚类
         type = "upper",
         col=colorRampPalette(c("green", "white", "red"))(50)
)
dev.off()

#第二个图
pdf(file="corpot2.pdf",width=8,height=8)
corrplot(M,
         order="original",
         method = "color",
         number.cex = 0.7, #相关系数
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

#第三个图
pdf(file="corpot3.pdf",width=7,height=7)
corrplot(M)
dev.off()

#第四个图
pdf(file="corpot4.pdf",width=8,height=8)
corrplot(M, order = "AOE", type = "upper", tl.pos = "d",tl.cex = 0.45)  #对称轴字体大小
corrplot(M, add = TRUE, type = "lower", method = "number", order = "AOE",
         diag = FALSE, tl.pos = "n", cl.pos = "n",number.cex=0.5 )     #相关性系数大小
dev.off()

#第五个图
pdf(file="corpot5.pdf",width=8,height=8)
chart.Correlation(M,method = "spearman")
dev.off()

#第六个图
pdf(file="corpot6.pdf",width=8,height=8)
corrplot(M, method="ellipse",p.mat = res$p, sig.level = 0.2,order = "AOE", type = "upper", tl.pos = "d",tl.cex = 0.45)
corrplot(M, add = TRUE, p.mat = res$p, sig.level = 0.2,type = "lower", method = "number", order = "AOE",
         diag = FALSE, tl.pos = "n", cl.pos = "n",number.cex=0.5)
dev.off()

#第七个图
pdf(file="corpot.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()

