
library(vioplot)       #引用包 
inputFile="CIBERSORT-Results.txt"     #免疫细胞浸润的结果文件
setwd("")     #设置工作目录

#读取免疫细胞浸润文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#对样品进行分组
con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]           #获取对照组的数据
treatData=rt[treat,]       #获取实验组的数据
conNum=nrow(conData)          #对照组样品数目
treatNum=nrow(treatData)      #实验组样品数目
rt=rbind(conData,treatData)

#绘制小提琴图
outTab=data.frame()
pdf(file="vioplot.pdf", width=13, height=8)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63), ylim=c(min(rt), max(rt)+0.05),
     xlab="", ylab="Fraction", main="", 
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞进行循环，绘制小提琴图，对照组用蓝色表示，实验组用红色表示
for(i in 1:ncol(rt)){
	  if(sd(rt[1:conNum,i])==0){
	    rt[1,i]=0.00001
	  }
	  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
	    rt[(conNum+1),i]=0.00001
	  }
	  #提取对照组和实验组样品的数据
	  conData=rt[1:conNum,i]
	  treatData=rt[(conNum+1):(conNum+treatNum),i]
	  #绘制小提琴图
	  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
	  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  #统计检验
	  wilcoxTest=wilcox.test(conData,treatData)
	  p=wilcoxTest$p.value
	  if(p<0.05){
	      cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
		  outTab=rbind(outTab,cellPvalue)
	  }
	  mx=max(c(conData,treatData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  #在小提琴图的上方加上差异的pvalue
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}

#绘制图例
legend("topright", 
       c("Control", "Treat"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))
#在X轴加上免疫细胞的名称
text(seq(1,64,3),-0.04,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

#输出免疫细胞和p值表格文件
write.table(outTab,file="immuneDiff.txt",sep="\t",row.names=F,quote=F)


