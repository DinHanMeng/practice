

library(limma)      #���ð�
expFile="normalize.txt"     #���������ļ�
geneFile="gene.txt"         #�����б��ļ�
setwd("")     #���ù���Ŀ¼

#��ȡ���������ļ����������ݽ��д���
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#��ȡ�����б��ļ�, ��ȡ����������ı�����
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#�������������ı�����
out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="FRGexp.txt", sep="\t", quote=F, col.names=F)
