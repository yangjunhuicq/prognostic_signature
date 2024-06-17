args <- commandArgs(TRUE)
if(length(args)!=6)
{
        print("<TPMfile> <OS_result_file> <outdir> <outprefix> <Groupby> <TRUE/FALSE: filter pvalue_sig or not>")
        print("OS time and status should be included in TPMfile.")
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

print("read in OS_result_file")
OS_res=read.table(args[2],header=T,stringsAsFactors=F,check.names=F,sep="\t")
head(OS_res[,c(1:4)],n=3)
OS_res=OS_res[-1,]
#OS_res=na.omit(OS_res,cols=c("ClusterID","Groupby"))
head(OS_res[,c(1:4)],n=3)
if(args[6]){
 OS_res=OS_res[which(OS_res$Groupby==args[5] & OS_res$pvalue_sig=="yes"),]
}else{
 OS_res=OS_res[which(OS_res$Groupby==args[5]),]
}
head(OS_res[,c(1:4)],n=3)
sig_genes=unique(OS_res$ClusterID)
directions=c()
for(i in c(1:length(sig_genes)))
{
 if(OS_res[i,"HR1"]>1){directions=c(directions,-1)}#high expression has better survival, i.e. this gene is highly expressed in patients with longer OS
 else{directions=c(directions,1)}
}
print(t(rbind(sig_genes,directions)))
warnings()

print("assign a score 1 or -1")
#For each sample, assign a score 1 or -1 for each gene based on its relative expression (> or <= median) and whether the signature gene is associated with bad or good survival.
TPMout=matrix(nrow=1,ncol=ncol(TPMdata))
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
write.table(TPMdata_siggenes,file=paste0(args[4],".TPM.siggenes.xls"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(TPMout,file=paste0(args[4],".TPM.siggenes.score.xls"),row.names=T,col.names=T,quote=F,sep="\t")
dim(TPMout)
dim(data_for_suv[,c(1,2)])
all(colnames(TPMout)==rownames(data_for_suv))
all(sort(colnames(TPMout))==sort(rownames(data_for_suv)))
TPMout_add_OS=cbind(data_for_suv[,c(1,2)],t(TPMout))
head(TPMout_add_OS[,c(1:4)],n=3)

out=matrix(NA,ncol=26)
plots=list()
OS_ana <- function(siggenes, ori_cind=ori_cindex,plots=plots)
{
  print("OS_ana start")
  #print(siggenes)
  #print(ori_cind)
  #print(length(siggenes))
  #print(head(TPMout_add_OS[,siggenes],n=3))
  #print(class(TPMout_add_OS[,siggenes]))
  if(length(siggenes)==1)
  {
   TPMout_add_OS$score=as.numeric(TPMout_add_OS[,siggenes])
  }else{
   TPMout_add_OS$score=apply(TPMout_add_OS[,siggenes],1,function(x) {sum(as.numeric(x))})
  }
  #print("add column score")
  TPMout_add_OS$risk=rep("low_risk",nrow(TPMout_add_OS))
  med=median(TPMout_add_OS$score)
  TPMout_add_OS$risk[TPMout_add_OS$score<=med] <- "high_risk"
  if(length(unique(TPMout_add_OS$risk))==1){
    print("######Note: only one group, so change the criteria [<=median vs other] to [<median vs other]")
    TPMout_add_OS$risk=rep("low_risk",nrow(TPMout_add_OS))
    TPMout_add_OS$risk[TPMout_add_OS$score<med] <- "high_risk"
  }
  #print("add column risk")

  factor="risk"
  {
        tmp_mat=TPMout_add_OS[,c(1,2,which(colnames(TPMout_add_OS)==factor))]
        colnames(tmp_mat)=c("TIME","STATUS",factor)
        tmp_mat_use2=na.omit(tmp_mat,cols=c("TIME","STATUS"))
        print(paste0("number of patients with OS info: ",nrow(tmp_mat_use2)))

        tmp_mat_usef <<- tmp_mat_use2
        write.table(tmp_mat_usef,file=paste0(args[4],".temp_OS_input_info.xls"),row.names=T,col.names=T,quote=F,sep="\t",append=T)
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
        out =rbind(out,c(paste0(siggenes,collapse=":"),"risk","merged",pvalue,group[1],HR[1],CI.lower[1],CI.upper[1],rep("",8),c_index,c_index.ci_low,c_index.ci_hig,cindex2$c.index, cindex2$lower, cindex2$upper,psig,hrsig1,hrsig2,hrsig3))
        #print("out ok")
        if(cindex2$c.index>ori_cind)
        {
          gene_signature=siggenes
          new_cindex=cindex2$c.index
        }else{
          gene_signature=siggenes[c(1:(length(siggenes)-1))]
          new_cindex=ori_cind
        }
        #print(c(gene_signature,new_cindex))
        #print(out)

        print("survfit...median")
        my_sc <- survfit(Surv(tmp_mat_usef$TIME,tmp_mat_usef$STATUS)~ tmp_mat_usef[,3],data=tmp_mat_usef)
        print("my_sc succeed")#print(my_sc)
        legend_labs <- unlist(sort(unique(tmp_mat_usef[,3])))
        print("legend_labs succeed")#gg <- ggsurvplot(my_sc,main = "Survival curve", pval = TRUE, legend = "top", legend.title = factor, legend.labs = legend_labs, risk.table = TRUE, surv.median.line='v')
        gg <- ggsurvplot(my_sc, data=tmp_mat_usef, pval = TRUE, legend = "top", legend.title = factor, legend.labs = legend_labs, risk.table = TRUE, surv.median.line='v', pval.method=TRUE, conf.int=TRUE)
        print("ggsurvplot succeed")
        plots[[length(plots)+1]]=gg

        print(paste("plots length: ",length(plots),sep=""))

        return(list(gene_signature,new_cindex,out,plots))
  }
}

print("Generation of prognostic signature")

seed_gene=OS_res$ClusterID[which.max(OS_res$c_index_survcomp)] #get the genename with highest C-index. If there are ties, return the first one
max_cindex=max(OS_res$c_index_survcomp)
sig_genes_sorted=OS_res$ClusterID[sort(OS_res$c_index_survcomp,index.return=T,decreasing=T)$ix]
sig_genes_sorted=sig_genes_sorted[-1]

outdata=OS_ana(seed_gene,ori_cind=0,plots=plots)
seed_gene=outdata[[1]]
max_cindex=outdata[[2]]
out=outdata[[3]]
plots=outdata[[4]]
print(paste0("The current signature is ",paste0(seed_gene,collapse=":"),"; c_index is ",max_cindex))
colnames(out)=c('ClusterID',"Groupby","cohort",'logrank.pvalue','group1','HR1','HR1.lower95','HR1.upper95','group2','HR2','HR2.lower95','HR2.upper95','group3','HR3','HR3.lower95','HR3.upper95',"c_index_survival","c_index_survival.ci_low","c_index_survival.ci_hig","c_index_survcomp","c_index_survcomp.ci_low","c_index_survcomp.ci_hig","pvalue_sig","HR1_sig","HR2_sig","HR3_sig")
write.table(out,file=paste(args[4],".temp_OS_res.file.xls",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

out=matrix(NA,ncol=26)
for(g in sig_genes_sorted)
{
  print(paste0("The current signature is ",paste0(seed_gene,collapse=":"),"; c_index is ",max_cindex))
  siggenes=c(seed_gene,g)
  print(paste0("test adding ",g))
  outdata=OS_ana(siggenes,ori_cind=as.numeric(max_cindex),plots=plots)
  seed_gene=outdata[[1]]
  max_cindex=outdata[[2]]
  out=outdata[[3]]
  plots=outdata[[4]]
  print(paste("plots length: ",length(plots),sep=""))
  print(paste0("The current signature is ",paste0(seed_gene,collapse=":"),"; c_index is ",max_cindex))
}
warnings()
print("writing output")
warnings()
colnames(out)=c('ClusterID',"Groupby","cohort",'logrank.pvalue','group1','HR1','HR1.lower95','HR1.upper95','group2','HR2','HR2.lower95','HR2.upper95','group3','HR3','HR3.lower95','HR3.upper95',"c_index_survival","c_index_survival.ci_low","c_index_survival.ci_hig","c_index_survcomp","c_index_survcomp.ci_low","c_index_survcomp.ci_hig","pvalue_sig","HR1_sig","HR2_sig","HR3_sig")
write.table(out,file=paste(args[4],".temp_OS_res.file.xls",sep=""),col.names=T,row.names=F,sep="\t",quote=F,append=T)

num=length(plots)
x=floor(num/6)-1
y=num-(x+1)*6
pdf(paste(args[4],".temp_OS_res.file.pdf",sep=""),width=30,height=7)
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

