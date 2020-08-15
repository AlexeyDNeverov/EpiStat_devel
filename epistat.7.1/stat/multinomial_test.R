#!/usr/bin/env Rscript
## Load the EMT package:
#library(EMT)
library(XNomial)
args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
if(is.na(infile)){
	print("The input file with data in a two columns' table is required! The first column is a probabilities the second data counts.");
	quit();
}
tbl=read.table(infile,header=TRUE)
observed=tbl$X
prob=tbl$P
#out <- multinomial.test(observed, prob, MonteCarlo=TRUE, ntrial=1e+8)
#plotMultinom(out,showmax = 10000)

xmonte(observed, prob, detail=3,ntrials=1e+6)
