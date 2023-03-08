
#���ð�
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
pvalueFilter=0.05     #pֵ��������
qvalueFilter=1        #�������pֵ��������

#������ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("")       #���ù���Ŀ¼
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #��ȡ�����ļ�

#��ȡ��������, ����������ת��Ϊ����id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #ɾ������idΪNA�Ļ���
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg��������
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
#�������õĹ�������,������������Ľ��
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#������������Ľ��
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#������ʾͨ·����Ŀ
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#��״ͼ
pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#����ͼ
pdf(file="bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()
