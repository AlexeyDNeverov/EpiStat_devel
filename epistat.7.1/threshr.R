#!/usr/bin/env Rscript
library(threshr)
args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
if(is.na(infile)){
	print("The input file with data in a two columns' table is required! The first column is a probabilities the second data counts.");
	quit();
}
in_lines=readLines(infile,skipNul=TRUE)
u_vec=as.numeric(unlist(strsplit(in_lines[1], "\t")))
d_vec=as.numeric(unlist(strsplit(in_lines[2], "\t")))
cv=ithresh(data=d_vec, u_vec=u_vec)
#graph_fn=paste(infile,"ithresh.graph.eps",sep=".")
#postscript(graph_fn)
#plot(cv)
#dev.off()
#graph_fn=paste(infile,"ithresh.GPD.eps",sep=".")
#postscript(graph_fn)
#plot(cv, which_u = "best")
#dev.off()
cv_sum=summary(cv)
best_u_idx=cv_sum[1,5]
best_u_samples=cv$sim_vals[cv$sim_vals[,4]==best_u_idx,]
mean_qu=mean(best_u_samples[,1])
mean_beta=mean(best_u_samples[,2])
mean_ksi=mean(best_u_samples[,3])
print(paste("best_u_idx=",best_u_idx))
print(paste("best_u=",u_vec[best_u_idx]))
print(paste("mean_Qu=",mean_qu))
print(paste("mean_beta=",mean_beta))
print(paste("mean_ksi=",mean_ksi))