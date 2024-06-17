args <- commandArgs(TRUE)
if(length(args)!=4)
{
        print("<TPMfile> <DEfile> <outdir> <outprefix> <avg_log2FC> <p_val_adj>")
        print("A title line should be included in infile.")
        print("Author: Yang Junhui")
        q()
}


setwd(args[3])
library(survival)
library(survminer)
library(gridExtra)
library(survcomp)
if(length(grep('csv$',args[1]))!=0)
{
  data_for_suv=read.csv(args[1],header=T,row.names=1,stringsAsFactors=F)
}else{
  data_for_suv=read.table(args[1],header=T,row.names=1,stringsAsFactors=F,sep="\t")
}
data_for_suv2=data_for_suv

print(args)
#print(as.numeric(args[4]))
out=matrix(NA,ncol=26)
plots=list()
#options(tmp_mat_usef='tmp_mat_usef')
OS_ana <- function(siggenes=siggenes, cluster='All',plots=plots)
{
  print(siggenes)
  genes=intersect(siggenes,colnames(data_for_suv))
  print(paste("sig genes for ",siggenes," in TPM file wrong: ",genes,sep=""))
  siggenes_mod=sub("-",".",siggenes,perl=T)
  genes=intersect(siggenes_mod,colnames(data_for_suv))
  print(paste("sig genes for ",siggenes," in TPM file: ",genes,sep=""))
  if(length(genes)>0){
  if(length(genes)==1)
  {
   data_for_suv2$sig=data_for_suv[,genes]
  }else{
   data_for_suv2$sig=apply(data_for_suv[,genes],1,mean)
  }
  #data_for_suv3=data_for_suv2[which(data_for_suv2$ACTARM=="MPDL3280A"),]
  factor=genes
  {
	tmp_mat=data_for_suv2[,c(1,2,which(colnames(data_for_suv2)==factor))]
	colnames(tmp_mat)=c("TIME","STATUS",factor)
	tmp_mat_use2=na.omit(tmp_mat,cols=c("TIME","STATUS"))
	print(paste0("number of patients with OS info: ",nrow(tmp_mat_use2)))
	tmp_mat_use2$median=rep("high",nrow(tmp_mat_use2))
	tmp_mat_use2$tertile=rep("lowT12",nrow(tmp_mat_use2))
	tmp_mat_use2$quartile=rep("lowT1",nrow(tmp_mat_use2))
	med=median(tmp_mat_use2[,3])
	tmp_mat_use2$median[tmp_mat_use2[,3]<=med] <- "low"
	x=floor(nrow(tmp_mat_use2)/3)
	order=sort(tmp_mat_use2[,3],decreasing=T,index.return=T)$ix
	tmp_mat_use2$tertile[order[c(1:x)]] <- "highT3"
	x=floor(nrow(tmp_mat_use2)/4)
	order=sort(tmp_mat_use2[,3],decreasing=T,index.return=T)$ix
	tmp_mat_use2$quartile[order[c(1:x)]] <- "highT4"
	tmp_mat_use2$quartile[order[c((x+1):(2*x))]] <- "lowT3"
	tmp_mat_use2$quartile[order[c((2*x+1):(3*x))]] <- "lowT2"
	#tmp_mat_use2$TPM6=rep("low6",nrow(tmp_mat_use2))
	#tmp_mat_use2$TPM6[tmp_mat_use2[,3]>6] <- "high6"

	tmp_mat_usef <<- tmp_mat_use2
	print("coxph...")
	colnames(tmp_mat_usef)=c("TIME","STATUS",factor,"Median","Tertile","Quartile")
	fit <- coxph(Surv(TIME,STATUS) ~ Median,data=tmp_mat_usef)
	test=summary(fit)
	pvalue <- test$sctest[["pvalue"]]
	group <- rownames(test$conf.int)
	HR <- as.numeric(test$conf.int[,1])
	CI.lower <- as.numeric(test$conf.int[,3])
	CI.upper <- as.numeric(test$conf.int[,4])
	#C-index by survival. Concordance = (#all concordant pairs + #tied pairs/2)/(#total pairs including ties). 考虑了tied risk。
	c_index <- test$concordance[[1]]
	# C-index的95%置信区间
	c_index.ci_low = c_index - test$concordance[[2]]*1.96
	c_index.ci_hig = c_index + test$concordance[[2]]*1.96
	#C-index by survcomp. Concordance = #all concordant pairs/#total pairs ignoring ties. 忽略tied risk （在处理tied risk上，也就是当两个观测拥有相同的生存时间和相同的自变量X时，忽略tied risk）
	cindex2 <- concordance.index(predict(fit), surv.time=tmp_mat_usef$TIME, surv.event=tmp_mat_usef$STATUS, method="noether",na.rm=TRUE)
	psig="no"
	if(pvalue<0.05){psig="yes"}
	hrsig1="no"
	if((HR[1]<1 & CI.upper[1]<1) | (HR[1]>1 & CI.lower[1]>1)){hrsig1="yes"}
	hrsig2=""
	hrsig3=""
	out =rbind(out,c(genes,"Median","merged",pvalue,group[1],HR[1],CI.lower[1],CI.upper[1],rep("",8),c_index,c_index.ci_low,c_index.ci_hig,cindex2$c.index, cindex2$lower, cindex2$upper,psig,hrsig1,hrsig2,hrsig3))

	fit <- coxph(Surv(TIME,STATUS) ~ Tertile,data=tmp_mat_usef)
	test=summary(fit)
	pvalue <- test$sctest[["pvalue"]]
	group <- rownames(test$conf.int)
	HR <- as.numeric(test$conf.int[,1])
	CI.lower <- as.numeric(test$conf.int[,3])
	CI.upper <- as.numeric(test$conf.int[,4])
	c_index <- test$concordance[[1]]
	# C-index的95%置信区间
	c_index.ci_low = c_index - test$concordance[[2]]*1.96
	c_index.ci_hig = c_index + test$concordance[[2]]*1.96
	#C-index by survcomp. Concordance = #all concordant pairs/#total pairs ignoring ties. 忽略tied risk （在处理tied risk上，也就是当两个观测拥有相同的生存时间和相同的自变量X时，忽略tied risk）
	cindex2 <- concordance.index(predict(fit), surv.time=tmp_mat_usef$TIME, surv.event=tmp_mat_usef$STATUS, method="noether",na.rm=TRUE)
	psig="no"
	if(pvalue<0.05){psig="yes"}
	hrsig1="no"
	if((HR[1]<1 & CI.upper[1]<1) | (HR[1]>1 & CI.lower[1]>1)){hrsig1="yes"}
	hrsig2=""
	hrsig3=""
	out =rbind(out,c(genes,"Tertile","merged",pvalue,group[1],HR[1],CI.lower[1],CI.upper[1],rep("",8),c_index,c_index.ci_low,c_index.ci_hig,cindex2$c.index, cindex2$lower, cindex2$upper,psig,hrsig1,hrsig2,hrsig3))

	fit <- coxph(Surv(TIME,STATUS) ~ Quartile,data=tmp_mat_usef)
	test=summary(fit)
	pvalue <- test$sctest[["pvalue"]]
	group <- rownames(test$conf.int)
	HR <- as.numeric(test$conf.int[,1])
	CI.lower <- as.numeric(test$conf.int[,3])
	CI.upper <- as.numeric(test$conf.int[,4])
	c_index <- test$concordance[[1]]
	# C-index的95%置信区间
	c_index.ci_low = c_index - test$concordance[[2]]*1.96
	c_index.ci_hig = c_index + test$concordance[[2]]*1.96
	#C-index by survcomp. Concordance = #all concordant pairs/#total pairs ignoring ties. 忽略tied risk （在处理tied risk上，也就是当两个观测拥有相同的生存时间和相同的自变量X时，忽略tied risk）
	cindex2 <- concordance.index(predict(fit), surv.time=tmp_mat_usef$TIME, surv.event=tmp_mat_usef$STATUS, method="noether",na.rm=TRUE)
	psig="no"
	if(pvalue<0.05){psig="yes"}
	hrsig1="no"
	if((HR[1]<1 & CI.upper[1]<1) | (HR[1]>1 & CI.lower[1]>1)){hrsig1="yes"}
	hrsig2="no"
	if((HR[2]<1 & CI.upper[2]<1) | (HR[2]>1 & CI.lower[2]>1)){hrsig2="yes"}
	hrsig3="no"
	if((HR[3]<1 & CI.upper[3]<1) | (HR[3]>1 & CI.lower[3]>1)){hrsig3="yes"}
	out =rbind(out,c(genes,"Quartile","merged",pvalue,group[1],HR[1],CI.lower[1],CI.upper[1],group[2],HR[2],CI.lower[2],CI.upper[2],group[3],HR[3],CI.lower[3],CI.upper[3],c_index,c_index.ci_low,c_index.ci_hig,cindex2$c.index, cindex2$lower, cindex2$upper,psig,hrsig1,hrsig2,hrsig3))
	print(dim(tmp_mat_usef))
	print("survfit...median")
	my_sc <- survfit(Surv(tmp_mat_usef$TIME,tmp_mat_usef$STATUS)~ tmp_mat_usef[,4],data=tmp_mat_usef)
	print("my_sc succeed")#print(my_sc)
	legend_labs <- unlist(sort(unique(tmp_mat_usef[,4])))
	print("legend_labs succeed")#gg <- ggsurvplot(my_sc,main = "Survival curve", pval = TRUE, legend = "top", legend.title = factor, legend.labs = legend_labs, risk.table = TRUE, surv.median.line='v')
	gg <- ggsurvplot(my_sc, data=tmp_mat_usef, pval = TRUE, legend = "top", legend.title = factor, legend.labs = legend_labs, risk.table = TRUE, surv.median.line='v')
	print("ggsurvplot succeed")
	plots[[length(plots)+1]]=gg
	print("survfit...tertile")
	my_sc2 <- survfit(Surv(tmp_mat_usef$TIME,tmp_mat_usef$STATUS)~ tmp_mat_usef[,5],data=tmp_mat_usef)
	#print(my_sc2)
	legend_labs <- unlist(sort(unique(tmp_mat_usef[,5])))
	gg <- ggsurvplot(my_sc2, data=tmp_mat_usef, pval = TRUE, legend = "top", legend.title = factor, legend.labs = legend_labs, risk.table = TRUE, surv.median.line='v')#, pval.method = T #show the method name used for calculating the pvalue, that corresponds to survival curves' comparison (Log-rank)
	plots[[length(plots)+1]]=gg
	print("survfit...quartile")
	my_sc3 <- survfit(Surv(tmp_mat_usef$TIME,tmp_mat_usef$STATUS)~ tmp_mat_usef[,6],data=tmp_mat_usef)
	#print(my_sc3)
	legend_labs <- unlist(sort(unique(tmp_mat_usef[,6])))
	gg <- ggsurvplot(my_sc3, data=tmp_mat_usef, pval = TRUE, legend = "top", legend.title = factor, legend.labs = legend_labs, risk.table = TRUE, surv.median.line='v')#, pval.method = T #show the method name used for calculating the pvalue, that corresponds to survival curves' comparison (Log-rank)
	plots[[length(plots)+1]]=gg
	print(paste("plots length: ",length(plots),sep=""))
	#write.table(tmp_mat_usef,file=paste(args[4],".tmp.",cluster,".tmp.xls",sep=""),col.names=T,row.names=T,quote=F,sep="\t")
  }
  }else{print(paste("no genes for analysis: ",genes))}
	return(list(plots,out))
}



if(length(grep('csv$',args[2]))!=0)
{
  infile=read.csv(args[2],header=F,row.names=1,stringsAsFactors=F)
}else{
  infile=read.table(args[2],header=F,row.names=1,stringsAsFactors=F,sep="\t")
}

for(g in rownames(infile))
{
  siggenes=g
  outdata=OS_ana(siggenes,plots=plots)
  plots=outdata[[1]]
  out=outdata[[2]]
}
print(paste("plots length: ",length(plots),sep=""))
print(paste("out dim: ",dim(out),sep=""))
print("writing output")
colnames(out)=c('ClusterID',"Groupby","cohort",'logrank.pvalue','group1','HR1','HR1.lower95','HR1.upper95','group2','HR2','HR2.lower95','HR2.upper95','group3','HR3','HR3.lower95','HR3.upper95',"c_index_survival","c_index_survival.ci_low","c_index_survival.ci_hig","c_index_survcomp","c_index_survcomp.ci_low","c_index_survcomp.ci_hig","pvalue_sig","HR1_sig","HR2_sig","HR3_sig")
write.table(out,file=paste(args[4],".xls",sep=""),col.names=T,row.names=F,sep="\t",quote=F)


num=length(plots)
x=floor(num/6)-1
y=num-(x+1)*6
pdf(paste(args[4],".pdf",sep=""),width=30,height=7)
i=-1
if(x>0)
{
 for(i in c(0:x))
 {
 	arrange_ggsurvplots(list(plots[[i*6+1]],plots[[i*6+2]],plots[[i*6+3]],plots[[i*6+4]],plots[[i*6+5]],plots[[i*6+6]]),print = TRUE,ncol = 6,nrow = 1, risk.table.height=0.35)
 }
}
re_list=list()
if(y>0)
{
 for(j in c(1:y))
 {
  re_list[[length(re_list)+1]]=plots[[i*6+6+j]]
 }
 arrange_ggsurvplots(re_list,print = TRUE,ncol = 6,nrow = 1, risk.table.height=0.35)
}
if(num==6)
{
 i=0
 arrange_ggsurvplots(list(plots[[i*6+1]],plots[[i*6+2]],plots[[i*6+3]],plots[[i*6+4]],plots[[i*6+5]],plots[[i*6+6]]),print = TRUE,ncol = 6,nrow = 1, risk.table.height=0.35)
}
dev.off()






