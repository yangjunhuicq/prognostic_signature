args <- commandArgs(TRUE)
if(length(args)!=4)
{
        print("<TPMfile> <gene_signature_file> <outdir> <outprefix>")
        print("OS time and status should be included in TPMfile.")
        print("gene_signature_file: containing genes for prognosis. with title; columns are: genename, direction(1/-1). 1 means that gene is highly expressed in patients with longer OS")
        print("Author: Yang Junhui")
        q()
}

library(survival)
library(survminer)
library(gridExtra)
library(survcomp)

setwd(args[3])
print("read in TPMfile")
if(length(grep('csv$',args[1]))!=0)
{
  data_for_suv=read.csv(args[1],header=T,row.names=1,stringsAsFactors=F)
}else{
  data_for_suv=read.table(args[1],header=T,row.names=1,stringsAsFactors=F,sep="\t") #sample in row. Colnames are OS_time, OS_status, genename1, genename2, ...
}
data_for_suv=na.omit(data_for_suv,cols=c(1,2)) #remove patients without survival information
TPMdata=t(data_for_suv[,-c(1,2)]) #gene in row, sample in col
head(data_for_suv[,c(1:4)],n=3)
head(TPMdata[,c(1:4)],n=3)

print("read in gene_signature_file")
OS_res=read.table(args[2],header=T,stringsAsFactors=F,check.names=F,sep="\t")
sig_genes=OS_res[,1]
directions=OS_res[,2] #high expression has better survival, i.e. this gene is highly expressed in patients with longer OS

print("assign a score 1 or -1")
#For each sample, assign a score 1 or -1 for each gene based on its relative expression (> or <= median) and whether the signature gene is associated with bad or good survival.
TPMout=matrix(nrow=1,ncol=ncol(TPMdata))
print("ERROR: the following signature genes are not in TPMfile. You should remove them from the gene_signature_file.")
setdiff(sig_genes,intersect(sig_genes,rownames(TPMdata)))
TPMdata_siggenes=TPMdata[intersect(sig_genes,rownames(TPMdata)),]
head(TPMdata_siggenes[,c(1:4)],n=3)
dim(TPMdata_siggenes)
for(i in c(1:length(sig_genes)))
{
 g=sig_genes[i]
 #print(paste0(i,":",g))
 thisline=rep("0",ncol(TPMout))
 thisline[which(TPMdata[g,] <= median(TPMdata[g,]))]=1*directions[i]
 thisline[which(TPMdata[g,] > median(TPMdata[g,]))]=-1*directions[i]
 TPMout=rbind(TPMout,thisline)
}
colnames(TPMout)=colnames(TPMdata)
TPMout=TPMout[-1,]
rownames(TPMout)=sig_genes
head(TPMout[,c(1:4)],n=3)
write.table(TPMdata_siggenes,file=paste0(args[4],".TPM.siggenes.validation.xls"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(TPMout,file=paste0(args[4],".TPM.siggenes.score.validation.xls"),row.names=T,col.names=T,quote=F,sep="\t")
dim(TPMout)
dim(data_for_suv[,c(1,2)])
all(colnames(TPMout)==rownames(data_for_suv))
all(sort(colnames(TPMout))==sort(rownames(data_for_suv)))
TPMout_add_OS=cbind(data_for_suv[,c(1,2)],t(TPMout))
head(TPMout_add_OS[,c(1:4)],n=3)

out=matrix(NA,ncol=26)
plots=list()
OS_ana <- function(siggenes, TPMinfo=TPMinfo,plots=plots,reptime=1,perc='1')
{
  print("OS_ana start")
  #print(siggenes)
  #print(ori_cind)
  #print(length(siggenes))
  #print(head(TPMout_add_OS[,siggenes],n=3))
  #print(class(TPMout_add_OS[,siggenes]))
  if(length(siggenes)==1)
  {
   TPMinfo$score=TPMinfo[,siggenes]
  }else{
   TPMinfo$score=apply(TPMinfo[,siggenes],1,function(x) {sum(as.numeric(x))})
  }
  #print("add column score")
  TPMinfo$risk=rep("low_risk",nrow(TPMinfo))
  med=median(TPMinfo$score)
  TPMinfo$risk[TPMinfo$score<=med] <- "high_risk"
  if(length(unique(TPMinfo$risk))==1){
    print("######Note: only one group, so change the criteria [<=median vs other] to [<median vs other]")
    TPMinfo$risk=rep("low_risk",nrow(TPMinfo))
    TPMinfo$risk[TPMinfo$score<med] <- "high_risk"
  }
  #print("add column risk")

  factor="risk"
  {
        tmp_mat=TPMinfo[,c(1,2,which(colnames(TPMinfo)==factor))]
        colnames(tmp_mat)=c("TIME","STATUS",factor)
        tmp_mat_use2=na.omit(tmp_mat,cols=c("TIME","STATUS"))
        print(paste0("number of patients with OS info: ",nrow(tmp_mat_use2)))

        tmp_mat_usef <<- tmp_mat_use2
        write.table(tmp_mat_usef,file=paste0(args[4],".temp_OS_input_info.validation.xls"),row.names=T,col.names=T,quote=F,sep="\t",append=T)
        print("coxph...")
        fit <- coxph(Surv(TIME,STATUS) ~ risk,data=tmp_mat_usef)
        test=summary(fit)
        print("coxph ok")
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
        #print("out")
        out =rbind(out,c(paste0(siggenes,collapse=":"),"risk",paste(reptime,perc,sep="_"),pvalue,group[1],HR[1],CI.lower[1],CI.upper[1],rep("",8),c_index,c_index.ci_low,c_index.ci_hig,cindex2$c.index, cindex2$lower, cindex2$upper,psig,hrsig1,hrsig2,hrsig3))
        #print("out ok")

        print("survfit...median")
        my_sc <- survfit(Surv(tmp_mat_usef$TIME,tmp_mat_usef$STATUS)~ tmp_mat_usef[,3],data=tmp_mat_usef)
        print("my_sc succeed")#print(my_sc)
        legend_labs <- unlist(sort(unique(tmp_mat_usef[,3])))
        print("legend_labs succeed")#gg <- ggsurvplot(my_sc,main = "Survival curve", pval = TRUE, legend = "top", legend.title = factor, legend.labs = legend_labs, risk.table = TRUE, surv.median.line='v')
        gg <- ggsurvplot(my_sc, data=tmp_mat_usef, pval = TRUE, legend = "top", legend.title = factor, legend.labs = legend_labs, risk.table = TRUE, surv.median.line='v', pval.method=TRUE, conf.int=TRUE)
        print("ggsurvplot succeed")
        plots[[length(plots)+1]]=gg

        print(paste("plots length: ",length(plots),sep=""))

        return(list(out,plots))
  }
}

print("Validation of prognostic signature")

out=matrix(NA,ncol=26)
for(i in c(1:100))
{
  sample_cols=sample(c(1:ncol(TPMdata)),size=round(ncol(TPMdata)*0.6))
  samples_per60=TPMout_add_OS[sample_cols,]
  samples_per40=TPMout_add_OS[-sample_cols,]
  outdata=OS_ana(sig_genes,TPMinfo=samples_per60,plots=plots,reptime=i,perc='0.6')
  out=outdata[[1]]
  plots=outdata[[2]]
  print(paste("plots length: ",length(plots),sep=""))

  outdata=OS_ana(sig_genes,TPMinfo=samples_per40,plots=plots,reptime=i,perc='0.4')
  out=outdata[[1]]
  plots=outdata[[2]]
  print(paste("plots length: ",length(plots),sep=""))
}
warnings()
print("writing output")
colnames(out)=c('ClusterID',"Groupby","cohort",'logrank.pvalue','group1','HR1','HR1.lower95','HR1.upper95','group2','HR2','HR2.lower95','HR2.upper95','group3','HR3','HR3.lower95','HR3.upper95',"c_index_survival","c_index_survival.ci_low","c_index_survival.ci_hig","c_index_survcomp","c_index_survcomp.ci_low","c_index_survcomp.ci_hig","pvalue_sig","HR1_sig","HR2_sig","HR3_sig")
write.table(out,file=paste(args[4],".temp_OS_res.file.validation.xls",sep=""),col.names=T,row.names=F,sep="\t",quote=F,append=F)

num=length(plots)
x=floor(num/6)-1
y=num-(x+1)*6
pdf(paste(args[4],".temp_OS_res.file.validation.pdf",sep=""),width=30,height=7)
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
warnings()
