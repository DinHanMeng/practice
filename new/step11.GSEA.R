


#���ð�
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
R.utils::setOption( "clusterProfiler.download.method",'auto' )

expFile="normalize.txt"     #���������ļ�
gene="TNFSF10"                 #��������
gmtFile="c2.cp.kegg.symbols.gmt"     #�����ļ�
setwd("")      #���ù���Ŀ¼

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ�����������Ʒ
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#����Ŀ��������������Ʒ���з��飬�õ��ߵͱ������logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #�ͱ����������
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #�߱����������
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#����logFC�Ի����������
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#��ȡ�����ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)
	
#����GSEA������ͼ��
termNum=6     #չʾͨ·����Ŀ
if(nrow(kkTab)>=termNum){
	showTerm=row.names(kkTab)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
	pdf(file="GSEA.pdf", width=7.5, height=5.5)
	print(gseaplot)
	dev.off()
}
