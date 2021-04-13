

setwd("C:\\Users\\lexb4\\Desktop\\immuneGene\\24.immuneCor")                     #设置工作目录

#读取文件，并对样品取交集
TIMER = read.table("immuneEstimation.txt", row.names=1 ,header=T,sep="\t",check.names=F)   #读取TIMER免疫细胞文件
immuneScore = read.table("risk.txt", row.names=1 ,header=T,sep="\t",check.names=F)   #读取生存相关的immuneScore文件
immuneScore=t(immuneScore[,3:(ncol(immuneScore)-1)])
TIMER=t(TIMER)
group=sapply(strsplit(colnames(TIMER),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
TIMER=TIMER[,group==0]
colnames(TIMER)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(TIMER))
sameSample=intersect(colnames(TIMER),colnames(immuneScore))
TIMER1=TIMER[,sameSample]
immuneScore1=immuneScore[,sameSample]

#相关性检验
outTab=data.frame()
i="riskScore"
for(j in row.names(TIMER1)){
	     x=as.numeric(immuneScore1[i,])
	     y=as.numeric(TIMER1[j,])
	     corT=cor.test(x,y)
	     cor=sprintf("%.03f",corT$estimate)
	     pvalue=corT$p.value
	     z=lm(y~x)
       pval=0
			 if(pvalue<0.001){
				   pval=signif(pvalue,4)
				   pval=format(pval, scientific = TRUE)
			 }else{
				   pval=sprintf("%.03f",pvalue)
			 }
			 pdf(file=paste0(j,".pdf"),height=4.5,width=4.5)
			 plot(x,y, type="p",pch=16,col="blue",main=paste("Cor=",cor," (p=",pval,")",sep=""),
				    cex=0.8, cex.lab=1.2, cex.main=1.2,cex.axis=1,     #由左到右一次是点、坐标轴、标题、刻度字体大小
				    xlab="Risk score",xlim=c(0,10),
				    ylab=j)
			 lines(x,fitted(z),col=2)
			 dev.off()
}

