#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
if(is.na(infile)){
	print("The input file with table having two rows is required! The first row is reference counts the second is data counts.");
	quit();
}
tbl=read.table(infile,header=FALSE)
#print(tbl)
counts=matrix(as.numeric(unlist(tbl)),nrow=nrow(tbl))
print(counts)
#fisher.test(counts,hybrid = TRUE,simulate.p.value = TRUE, B = 2000000)
if(ncol(tbl)==2&&nrow(tbl)==2){
	x=try(fisher.test(counts,workspace = 1e+6))
}else{
	x=try(fisher.test(counts,workspace = 1e+6,hybrid = TRUE))
}
if(inherits(x, "try-error")){
	fisher.test(counts,hybrid = TRUE,simulate.p.value = TRUE, B = 1e+8)
}else{
	print(x)
}
