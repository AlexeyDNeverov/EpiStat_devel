#!/usr/bin/env Rscript
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite(c("BiRewire"))
#.libPaths(c( .libPaths(), "c:/Users/neva/Documents/R/win-library/3.2"))
library(BiRewire)
args <- commandArgs(trailingOnly = TRUE)
inMatrixFN=args[1]
outFolder=args[2]
if(is.na(outFolder)){
	outFolder="./"
}else{
	if(regexpr("/$",outFolder)==-1){
		outFolder=paste(outFolder,"/", sep="")
	}
}
from=as.numeric(args[3])
if(is.na(from)){
	from=1
}
to=as.numeric(args[4])
if(is.na(to)||to<from){
	to=from
}
#input matrix: branches in rows; sites in columns
#output: rows - list of sites carrying mutations on a particular branch
data <- read.table(inMatrixFN, sep="", colClasses='character')
dat <-sapply(data, function(e) {
  splitter <- strsplit(e,"")
  
})
dat <- sapply(dat, function (e) {
  as.numeric(unlist(e))
})
dat <- t(dat)
for (i in from:to){
	m2<-birewire.rewire.bipartite(dat,verbose=FALSE)
	subs_on_node <-apply(m2, 1, function(e){
		which(e>0)
	})
	tv<-sapply(subs_on_node, function(e){
		paste(e, collapse=";")
	})
	cat (tv, sep="\n",file=paste(outFolder,i,".mlist", sep=""))
}