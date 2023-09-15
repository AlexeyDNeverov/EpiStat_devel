#!/usr/bin/env Rscript
library("SuperExactTest")
args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
outfile_pref=args[2]
if(is.na(infile)){
	print("The input file with data is required!");
	quit();
}
if(is.na(outfile_pref)){
	outfile_pref=infile;
}
in_lines=readLines(infile,skipNul=TRUE)
n_total=as.numeric(in_lines[1])
n_sets=length(in_lines)-1
x=vector("list",n_sets)
for(i in 1:n_sets){
	x[[i]]=as.numeric(unlist(strsplit(in_lines[i+1], "\t")))
}
res=supertest(x, n=n_total)
graph_fn=paste(outfile_pref,"SuperExactTest.eps",sep=".")
postscript(graph_fn)
plot(res, sort.by="size")
dev.off()
summary(res)
