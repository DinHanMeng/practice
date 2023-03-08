

library(limma)      #引用包
expFile="normalize.txt"     #表达数据文件
geneFile="gene.txt"         #基因列表文件
setwd("")     #设置工作目录

#读取表达数据文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取铁死亡基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#输出铁死亡基因的表达量
out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="FRGexp.txt", sep="\t", quote=F, col.names=F)

