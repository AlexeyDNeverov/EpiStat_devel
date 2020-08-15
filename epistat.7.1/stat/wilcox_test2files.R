#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
infile1=args[1]
infile2=args[2]
altern_opt=args[3]
if(is.na(altern_opt)){
	altern_opt="two.sided"
}
tbl1=read.delim(infile1)
x=na.omit(tbl1[[ncol(tbl1)]])
cat("nx=",length(x),"Mx=",mean(x),"Sx=",var(x)**0.5,"\n") 
tbl2=read.delim(infile2)
y=na.omit(tbl2[[ncol(tbl2)]])
cat("ny=",length(y),"My=",mean(y),"Sy=",var(y)**0.5)
wilcox.test(x,y,alternative=altern_opt)