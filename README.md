1. step1.survival.notitle.cindex.R
1) arguments:
	<TPMfile>	input file. The first column is samplename, the second column is survival time, the third column is survival status. From the 4th column are genesymbols.
	<genelist>	input file without header. One gene in each row.
	<outdir>	output directory.
	<outprefix>	output prefix.
2) output files:
	outprefix.xls	result file for each gene.
	outprefix.pdf	KM plot

2. step2.prognostic_signature.and.validation.v4.sh
1) arguments:
	<TPMfile>	input file. The first column is samplename, the second column is survival time, the third column is survival status. From the 4th column are genesymbols.
	<OS_result_file>	outprefix.xls from step1.
	<outdir>	output directory.
	<outprefix>	output prefix.
	<Groupby>	Median/Tertile/Quartile
	<TRUE/FALSE: filter pvalue_sig or not>
	<infile>  TPMfile for validation cohort
2) output files:
	outprefix.temp_OS_input_info.xls  result file for each gene signature.
	outprefix.temp_OS_res.file.pdf  KM plot.
	outprefix.gene_signature.xls  The prognostic signature. One gene per row.
	outprefix.validation.temp_OS_res.file.validation.xls  The OS result for 100 validations. 
	outprefix.validation.temp_OS_res.file.validation.pdf  KM plot for 100 validations. 