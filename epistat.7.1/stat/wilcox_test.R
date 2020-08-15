#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
altern_opt=args[2]
if(is.na(altern_opt)){
	altern_opt="two.sided"
}
tbl=read.delim(infile)
header=colnames(tbl)
n=ncol(tbl)
f=paste(header[c(n-1,n)],collapse=' ~ ')
wilcox.test(formula=as.formula(f),data=tbl,alternative=altern_opt)