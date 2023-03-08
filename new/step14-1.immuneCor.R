
#引用包
library(limma)
library(reshape2)
library(ggpubr)
library(ggExtra)

gene="MAPK8"                       #基因名称
expFile="normalize.txt"             #表达数据文件
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
setwd("")     #设置工作目录

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#获取目标基因表达量
data=t(data[gene,,drop=F])
data=as.data.frame(data)

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])

#对免疫细胞进行循环，绘制相关性散点图
outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-1)]){
	x=as.numeric(rt[,gene])
	y=as.numeric(rt[,i])
	if(sd(y)==0){y[1]=0.00001}
	cor=cor.test(x, y, method="spearma")
	
	outVector=cbind(Gene=gene, Cell=i, cor=cor$estimate, pvalue=cor$p.value)
	outTab=rbind(outTab,outVector)
	
	if(cor$p.value<0.05){
		outFile=paste0("cor.", i, ".pdf")
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				  xlab(paste0(gene, " expression")) + ylab(i)+
				  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				  stat_cor(method = 'spearman', aes(x =x, y =y))
		p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
		#相关性图形
		pdf(file=outFile, width=5.2, height=5)
		print(p2)
		dev.off()
	}
}
#输出免疫功能和p值表格文件
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)


