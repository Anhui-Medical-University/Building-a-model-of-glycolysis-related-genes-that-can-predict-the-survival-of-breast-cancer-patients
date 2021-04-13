
#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

setwd("C:\\Users\\Administrator\\Desktop\\96glycolysis\\20.cliSurvival")                 #设置工作目录
risk=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)        #读取风险文件
cli=read.table("clinical.txt",sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
data=cbind(futime=risk[,1],fustat=risk[,2],cli)

for(i in colnames(data[,3:ncol(data)])){
    rt=data[,c("futime","fustat",i)]
    rt=rt[(rt[,i]!="unknow"),]
    colnames(rt)=c("futime","fustat","clinical")
	tab=table(rt[,"clinical"])
	tab=tab[tab!=0]
	labels=paste0(names(tab),"(n=",tab,")")
	diff=survdiff(Surv(futime, fustat) ~clinical,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ clinical, data = rt)
	#绘制生存曲线
	surPlot=ggsurvplot(fit, 
	           data=rt,
	           pval=pValue,
	           pval.size=6,
	           legend.labs=labels,
	           legend.title=i,
	           xlab="Time(years)",
	           break.time.by = 1,
	           palette=c("blue","red") )
    #输出图片
	pdf(file=paste0("survival.",i,".pdf"),onefile = FALSE,
	       width = 5,             #图片的宽度
	       height =4.5)           #图片的高度
	print(surPlot)
	dev.off()
}

