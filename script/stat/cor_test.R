#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
if(is.na(infile)){
	print("The input file with data in a two columns' table is required! The first column is a probabilities the second data counts.");
	quit();
}
method=args[2]
if(is.na(method)){
	warning("Error: correlation statistic is not specified! The Pearson's statistic is used.")
	method="pearson"
}
if(!(method=="spearman"||method=="pearson"||method=="kendall")){
	warning(paste("Error: unexpected correlation statistic:",method))
	quit()
}
tbl=read.table(infile,header=TRUE)
cor.test(tbl[[1]],tbl[[2]],method=method)