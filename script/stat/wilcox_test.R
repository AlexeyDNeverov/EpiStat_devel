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
x=subset(tbl,select=header[n-1],subset=(tbl[,n]==1))
y=subset(tbl,select=header[n-1],subset=(tbl[,n]!=1))
paste(c("#",header[n],"=",sum(!is.na(x[,1]))),collapse="")
paste(c("mean(",header[n],")=",mean(x[,1], na.rm=TRUE)),collapse="")
paste(c("#complement=",sum(!is.na(y[,1]))),collapse="")
paste(c("mean(complement)=",mean(y[,1], na.rm=TRUE)),collapse="")