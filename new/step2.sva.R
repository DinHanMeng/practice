
#���ð�
library(limma)
library(sva)
outFile="data.txt"       #����ļ�
setwd("")      #���ù���Ŀ¼

#��ȡĿ¼������".txt"��β���ļ�
files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#��ȡ����txt�ļ��еĻ�����Ϣ�����浽geneList
for(file in files){
	if(file==outFile){next}
    rt=read.table(file, header=T, sep="\t", check.names=F)      #��ȡ�����ļ�
    geneNames=as.vector(rt[,1])      #��ȡ��������
    uniqGene=unique(geneNames)       #����ȡunique
    header=unlist(strsplit(file, "\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}

#��ȡ��������
interGenes=Reduce(intersect, geneList)

#���ݺϲ�
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    #��ȡ�����ļ������������ļ���������
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)

    #����ֵ�������ȡlog2
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    rt=normalizeBetweenArrays(rt)
    
    #���ݺϲ�
    if(i==1){
    	allTab=rt[interGenes,]
    }else{
    	allTab=cbind(allTab, rt[interGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#�����ݽ��н��������������Ľ��
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)

