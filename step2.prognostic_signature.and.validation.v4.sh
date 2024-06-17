TPMfile=$1
OS_result_file=$2
outdir=$3
prefix=$4
Groupby=$5
filter_psig=$6 #with no quote
infile=$7 #TPMfile for validation cohort
finalprefix=$prefix.validation
Rscript prognostic_signature.v4.R $TPMfile $OS_result_file $outdir $prefix $Groupby $filter_psig >>$prefix.prognostic_signature.log 2>>$prefix.prognostic_signature.log
echo "gene	directions" >$prefix.gene_signature.xls
grep -v -E "NA|c_index_survcomp" $prefix.temp_OS_res.file.xls |awk -F"\t" 'BEGIN{max_cindex=0;gene_signature=""}{if(NR>1 && $20>max_cindex){max_cindex=$20;gene_signature=$1}}END{print gene_signature"\t"max_cindex}'|cut -f1|awk -F":" '{for(i=1;i<=NF;i++)print "\""$i"\""}'|grep -f - -w $prefix.prognostic_signature.log|sed 's/^ \+//;s/\[.\+\] \+//;s/ \+/	/;s/"//g' >>$prefix.gene_signature.xls
Rscript prognostic_signature.validate.v3.R $infile $prefix.gene_signature.xls $outdir $finalprefix >>$finalprefix.validation.log 2>>$finalprefix.validation.log
