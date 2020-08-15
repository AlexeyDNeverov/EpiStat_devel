#!/usr/bin/env Rscript
.libPaths(c( .libPaths(), "~/R/3.2/RPackages"))
#.libPaths(c( .libPaths(), "C:\\Users\\neva\\Documents\\R\\win-library\\3.2"))
library('corpcor')
library('getopt')

spec = matrix(c(
	'file','f',1,"character","(Required) Block matrix file name",
	'lambda','l',1,"double","Shrincage intensity [0,1] of the input correlation matrix. Default=0 (estimate)",
	'suffix','s',1,"character","Suffix of output file with inverse matrix",
	'norm','n',1,"integer","Normalization mode of inverse matrix:\n\t0 - no normalization,\n\t1 - average absolute value,\n\t2 (Default) - APC,\n\t3 - max absolute value",
	'transform','t',1,"integer","Transformation mode of an input correlation matrix:\n\t0 (Default) - no transformation,\n\t1 - absolute values,\n\t2 - geometric mean of max abs row and max abs col",
	'help','h',0,"logical","This help message",
	'keep_unchanged','k',0,"logical","Do not invert input correlation matrix, only normalization and transformations are applied"
	),byrow=TRUE, ncol=5);
opt=getopt(spec)
if(!is.null(opt$help)){
	cmd="COR2PCOR - calculates partial correlations for a specified correlation matrix: <block_matrix>\n"
	cmd=paste(cmd,get_Rscript_filename())
	cat(getopt(spec,command=cmd, usage=TRUE))
	q(status=1);
}
if(is.null(opt$file)){
	write("The name of an input file is required! Use --file <file_name> option",stderr())
	q(status=1)
}
if(is.null(opt$lambda)){opt$lambda=0}
if(is.null(opt$suffix)){opt$suffix=".cor2pcor.R.out"}
if(is.null(opt$norm)){opt$norm=2}
if(is.null(opt$transf)){opt$transform=0}
if(is.null(opt$keep_unchanged)){opt$keep_unchanged=FALSE}

inMatrixFN=opt$file
dens=opt$lambda
outFNExt=opt$suffix
message("input matrix: ", inMatrixFN)
	
readSubMatrix<-function(inMatrixFN){
	con  <- file(inMatrixFN, open = "r")
	y1=-1
	while(y1==-1) {
		if(length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0){
			y1=regexpr("Matrix:\\s*\\S*?\\[",oneLine,perl = TRUE)
		}
	} 
	x=read.csv(con,sep='\t')
	close(con)
	x
}

#Reading block matrix file
D=0
con  <- file(inMatrixFN, open = "r")
while(D==0){
	if(length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0){
		y1=regexpr("Block Matrix:\\s*\\S*?\\[",oneLine,perl = TRUE)
		if(y1!=-1){
			d1=0
			d2=0
			z1 <- y1 + attr(y1, "match.length")
			oneLine=substr(oneLine,z1,nchar(oneLine))
			y1=regexpr("\\d+",oneLine,perl = TRUE)
			if(y1!=-1){
				z1 <- y1 + attr(y1, "match.length")-1
				d1=as.numeric(substr(oneLine,y1,z1))
				oneLine=substr(oneLine,z1+1,nchar(oneLine))
			}
			y1=regexpr("\\d+",oneLine,perl = TRUE)
			if(y1!=-1){
				z1 <- y1 + attr(y1, "match.length")-1
				d2=as.numeric(substr(oneLine,y1,z1))
			}
			if(d1>0 & d2==d1){
				D=d1
			}else if(d1>0){
				message("The block matrix is not symmetric!")
			}else{
				message("Unable to identify matrix dimensions!")
			}
		}
	}
}
inputDir=""
submtxExt=""
if(D){
	repeat{
		if(length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0){
			y1=regexpr("Submatrices:",oneLine,perl = TRUE)
			if(y1!=-1) break
			y1=regexpr('OutDir\\s*=\\s*"',oneLine,perl = TRUE)
			if(y1!=-1){
				z1 <- y1 + attr(y1, "match.length")
				oneLine=substr(oneLine,z1,nchar(oneLine))
				y1=regexpr('"',oneLine,perl = TRUE)
				if(y1!=-1){
					z1 <- y1-1
					inputDir=substr(oneLine,1,z1)
				}
			}
			y1=regexpr('SubMtxExt\\s*=\\s*"',oneLine,perl = TRUE)
			if(y1!=-1){
				z1 <- y1 + attr(y1, "match.length")
				oneLine=substr(oneLine,z1,nchar(oneLine))
				y1=regexpr('"',oneLine,perl = TRUE)
				if(y1!=-1){
					z1 <- y1-1
					submtxExt=substr(oneLine,1,z1)
				}
			}
		}
	}
}
x=read.csv(con,sep='\t',na.strings=".",stringsAsFactors = FALSE)
close(con)
#initializing input matrix
sub_mtx_data_list=list()
if(D*D==nrow(x)){
	sub_mtx_names=vector()
	k=1
	for(i in 1:(D*D)){
		name=x[i,]$name
		if(name!="0" & name!="I"){
			sub_mtx_names[k]=name
			k=k+1
		}
	}
	sub_mtx_names=unique(sub_mtx_names)
	for(i in 1:length(sub_mtx_names)){
		name=sub_mtx_names[i]
		mtx_fn=paste(inputDir,name,submtxExt, sep = "", collapse = "")
		sub_mtx_data_list[[name]]=readSubMatrix(mtx_fn)
	}
}else{
	message("Wrong number of blocks!")
	D=0
}

diagMtx<-function(nRows,nCols){diag(1,nrow=nRows,ncol=nCols)}
zeroMtx<-function(nRows,nCols){matrix(0,nrow=nRows,ncol=nCols)}

mkSymmetric<-function(mtx){
	d=nrow(mtx)
	M=mtx
	for(i in 2:d){
		for(j in 1:(i-1)){
			if(abs(mtx[i,j])>abs(mtx[j,i])){
				M[j,i]=mtx[i,j]
			}else{
				M[i,j]=mtx[j,i]
			}
		}
	}
	M
}

absTransform<-function(mtx){
	d=nrow(mtx)
	#M=matrix(abs(mtx),nrow=d,ncol=d)
	M=abs(mtx)
}

geomeanMaxR_MaxC_Transform<-function(mtx){
	d=nrow(mtx)
	row_max=sapply(1:d, function(i, y, z){ A=0; B=0; if(i-1>0) A=max(abs(y[i,1:(i-1)])); if(i+1<=d) B=max(abs(y[i,(i+1):d])); max(A,B)}, y=mtx, z=d)
	col_max=sapply(1:d, function(i, y, z){ A=0; B=0; if(i-1>0) A=max(abs(y[1:(i-1),i])); if(i+1<=d) B=max(abs(y[(i+1):d,i])); max(A,B)}, y=mtx, z=d)
	M=mtx
	for(i in 1:d){
		if(row_max[i]==0){
			message("Error - empty row: ", i)
			quit()
		}
		for(j in 1:d){
			if(col_max[j]==0){
				message("Error - empty column: ", j)
				quit()
			}
			if(i==j){
				M[i,j]=1
			}else{
				M[i,j]=M[i,j]/(sqrt(row_max[i])*sqrt(col_max[j]))
			}
		}
	}
	M
}

shrinkMatrix<-function(mtx,lambda,nIter){
	dmean=0
	d=nrow(mtx)
	#dMtx=d/2
	dmean=0
	for(i in 1:d){
		dmean=dmean+mtx[i,i]
	}
	dmean=dmean/d
	dmean
	M=mtx
	x=try(chol(M))
	I=1
	while(inherits(x, "try-error") & I<=nIter){
		for(i in 1:d){
			for(j in 1:d){
				if(i==j){
					M[i,j]=dmean*lambda+(1.0-lambda)*M[i,j]
				}else{
					M[i,j]=(1.0-lambda)*M[i,j]
				}
			}
		}
		x=try(chol(M))
		print(I)
		I=I+1
	}
	M
}

estimateIvcov<-function(mtx,iDens){
	if(iDens>0){
		X=pcor.shrink(mtx,iDens)
	}else{
		X=pcor.shrink(mtx)
	}
	X
}

mkNormResults<-function(mtx,nmode){
	dMtx=nrow(mtx)
	row_margs=rep(0,dMtx)
	col_margs=rep(0,dMtx)
	mtx_margs=0
	mtx_max=0
	for(i in 2:dMtx){
        for(j in 1:(i-1)){
			if(i != j){
				if(!is.na(mtx[i,j])){
					pc=abs(mtx[i,j])
					if(mtx_max<pc){mtx_max=pc}
					row_margs[i]=row_margs[i]+pc
					col_margs[j]=col_margs[j]+pc
					mtx_margs=mtx_margs+pc
				}
			}
		}
	}
	mtx_margs=mtx_margs/( dMtx*(dMtx-1)/2 )
	message("mkNormResults() mtx_margs: ", mtx_margs)
	for(i in 1:dMtx){
		tmp = row_margs[i]+col_margs[i]
		col_margs[i]=tmp
		row_margs[i]=tmp
		row_margs[i]=row_margs[i]/(dMtx-1)
		col_margs[i]=col_margs[i]/(dMtx-1)
	}
	M=mtx
	for(i in 1:dMtx){
		for(j in 1:dMtx){
			if(is.na(mtx[i,j])){
				M[i,j]=0
			}else if(nmode==1){
				M[i,j]=M[i,j]/mtx_margs
			}else if(nmode==2){
				M[i,j]=(abs(M[i,j])-(row_margs[i]*col_margs[j])/mtx_margs)
			}else if(nmode==3){
				M[i,j]=M[i,j]/mtx_max
			}else if(nmode>2){
				write("Error in mkNormResults(): Unknown normalization mode",stderr())
				q(status=1)
			}
		}
	}
	M
}

mateMtx<-function(invMtx,mtx){
	dMtx=nrow(mtx)
	M=invMtx
	for(i in 2:dMtx){
        for(j in 1:(i-1)){
			w=invMtx[i,j]
			if(abs(invMtx[i,j])<abs(invMtx[j,i])){
				w=invMtx[j,i]
			}
			if(mtx[i,j]!=mtx[j,i]){
				if(abs(mtx[i,j])>abs(mtx[j,i])){
					M[i,j]=w
					M[j,i]=0
				}else{
					M[j,i]=w
					M[i,j]=0
				}
			}else if(mtx[i,j]==0){
				M[j,i]=0
				M[i,j]=0
			}
		}
	}
	M
}

#start script
if(D>0){
	mcol=list()
	for(i in 1:D){
		mrow=list()
		nRows=x[1+(i-1)*D,]$nrow
		for(j in 1:D){
			k=j+(i-1)*D
			name=x[k,]$name
			nCols=x[k,]$ncol
			if(name == "0"){
				C=zeroMtx(nRows,nCols)
			}else if(name == "I"){
				C=diagMtx(nRows,nCols)
			}else{
				C=matrix(sub_mtx_data_list[[name]][[3]],nRows,nCols,byrow = TRUE)
			}
			mrow[[j]]=C
		}
		C=mrow[[1]]
		j=2
		while(j<=D){
			if(nrow(C)==nrow(mrow[[j]])){
				C=cbind(C,mrow[[j]])
			}else{
				message("Wrong number of rows in the block number: ",j+(i-1)*D)
			}
			j=j+1
		}
		mcol[[i]]=C
	}
	C=mcol[[1]]
	i=2
	while(i<=D){
		if(ncol(C)==ncol(mcol[[i]])){
			C=rbind(C,mcol[[i]])
		}else{
			message("Wrong number of columns in row: ",i)
		}
		i=i+1
	}
	mcol=list()
	if(opt$transform==1){
		C=absTransform(C)
	}else if(opt$transform==2){
		C=geomeanMaxR_MaxC_Transform(C)
	}else if(opt$transform>2){
		write("Error: Unknown input matrix transformation mode",stderr())
		q(status=1)
	}
	X=mkSymmetric(C)
	if(opt$keep_unchanged){
		write("Warning message:\nThe input matrix wasn't inverted because the '--keep_unchanged (-k)' option accounted!",stderr())
	}else{
		X=estimateIvcov(X,dens)
	}
	wwi=mkNormResults(X,opt$norm)
	wwi=mateMtx(wwi,C)
	d1=0
	for(i in 1:D){
		nRows=x[1+(i-1)*D,]$nrow
		d1=d1+nRows
		d2=0
		for(j in 1:D){
			k=j+(i-1)*D
			name=x[k,]$name
			nCols=x[k,]$ncol
			d2=d2+nCols
			if(name!="0" & name!="I" & is.na(x[k,]$operator)){
				sub_mtx_data_list[[name]]$wi_score=as.vector(t(wwi[(d1-nRows+1):d1,(d2-nCols+1):d2]))
				out_fn=paste(inputDir,name,outFNExt, sep = "", collapse = "")
				write.table(sub_mtx_data_list[[name]],out_fn,quote=FALSE, row.names = F, sep='\t')
			}
		}
	}
}

