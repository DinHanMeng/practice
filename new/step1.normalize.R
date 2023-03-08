
library(limma)               #���ð�
expFile="merge.txt"     #���������ļ�
conFile="s1.txt"             #��������Ʒ�ļ�
treatFile="s2.txt"           #ʵ������Ʒ�ļ�
setwd("")      #���ù���Ŀ¼

#��ȡ�����ļ������������ļ�����
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#��ȡ��������Ʒ�ļ�,��ȡ��������Ʒ�ı�������
s1=read.table(conFile, header=F, sep="\t", check.names=F)
sampleName1=as.vector(s1[,1])
conData=data[,sampleName1]

#��ȡʵ������Ʒ�ļ�,��ȡʵ������Ʒ�ı�������
s2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName2=as.vector(s2[,1])
treatData=data[,sampleName2]

#���ݺϲ�
rt=cbind(conData, treatData)

#�������û��ȡlog2,��������Զ�ȡlog2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	rt[rt<0]=0
	rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

#���������ı�����, ͬʱ����Ʒ���ֺ��������Ʒ�ķ�����Ϣ
conNum=ncol(conData)
treatNum=ncol(treatData)
Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)