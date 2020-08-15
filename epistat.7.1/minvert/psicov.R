#!/usr/bin/env Rscript
.libPaths(c( .libPaths(), "~/R/3.2/RPackages"))
#.libPaths(c( .libPaths(), "C:\\Users\\neva\\Documents\\R\\win-library\\3.2"))
library('glasso')
library('getopt')

spec = matrix(c(
	'file','f',1,"character","(Required) Block matrix file name",
	'density','d',1,"double","Desired dencity of inverse matrix",
	'suffix','s',1,"character","Suffix of output file with inverse matrix",
	'norm','n',1,"integer","Normalization mode of inverse matrix:\n\t0 - no normalization,\n\t1 - average absolute value,\n\t2 (Default) - APC,\n\t3 - max absolute value",
	'transform','t',1,"integer","Transformation mode of an input correlation matrix:\n\t0 (Default) - no transformation,\n\t1 - absolute values",
	'help','h',0,"logical","This help message"
	),byrow=TRUE, ncol=5);
opt=getopt(spec)
if(!is.null(opt$help)){
	cmd="PSICOV - performes graphycal lasso for inversion of a specified correlation matrix: <block_matrix>\n"
	cmd=paste(cmd,get_Rscript_filename())
	cat(getopt(spec,command=cmd, usage=TRUE))
	q(status=1);
}
if(is.null(opt$file)){
	write("The name of an input file is required! Use --file <file_name> option",stderr())
	q(status=1)
}
if(is.null(opt$density)){opt$density=0}
if(is.null(opt$suffix)){opt$suffix=".psicov.R.out"}
if(is.null(opt$norm)){opt$norm=2}
if(is.null(opt$transf)){opt$transform=0}

inMatrixFN=opt$file
dens=opt$density
outFNExt=opt$suffix
message("mkNormResults() inMatrixFN: ", inMatrixFN)
	
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
message("mkNormResults() reading block matrix: ", inMatrixFN)
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
	#d=nrow(mtx)
	#M=matrix(abs(mtx),nrow=d,ncol=d)
	M=abs(mtx)
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

matrixDensity<-function(mtx){
    nz = 0
	dens=0
	d=nrow(mtx)
    for(i in 2:d){
        for(j in 1:(i-1)){
			if(mtx[i,j]!= 0){
				nz=nz+1
			}
		}
	}
    dens = nz / ( d * (d-1)/2 )
}

mkRhoMatrix<-function(dMtx,dis_idx,rho,BIG){
	rhoMtx=matrix(1,nrow=dMtx,ncol=dMtx)
	rhoMtx=rho*rhoMtx
	k=nrow(dis_idx)
	for(i in 1:k){
		rhoMtx[dis_idx[i,1],dis_idx[i,2]]=BIG
		rhoMtx[dis_idx[i,2],dis_idx[i,1]]=BIG
	}
	rhoMtx
}

estimateIvcov<-function(mtx,iDens){
	BIG=1000000
	rho = 0.001
	lambda=.2
	d=nrow(mtx)
	nIter=100
	X=shrinkMatrix(mtx,lambda,nIter)
	if(iDens==0){
		a=glasso(X,rho=0)
		wi=a$wi
		return(wi)
	}
	k=0
	for(i in 2:d){
        for(j in 1:(i-1)){
			if(mtx[i,j]==0){
				k=k+1
			}
		}
	}
	dis_idx=matrix(0,nrow=k,ncol=2)
	k=1
	for(i in 2:d){
        for(j in 1:(i-1)){
			if(mtx[i,j]==0){
				dis_idx[k,1]=i
				dis_idx[k,2]=j
				k=k+1
			}
		}
	}
	k=k-1
	message("estimateIvcov() rho: ", rho)
	rhoMtx=mkRhoMatrix(d,dis_idx,rho,BIG)
	a=glasso(X,rho=rhoMtx,zero=dis_idx)
#a=glasso(X,rho=rho)
	wi=a$wi
	prev_dens=-1
	dens=matrixDensity(wi)
	rhofac=0
	denseDiff=abs(iDens - dens)/dens
	bestDenseDiff=denseDiff
	bestRho=rho
	bestDense=dens
	I=1
	while(I<=nIter & (dens == 0 | denseDiff > 0.01) ){
		message("estimateIvcov() matrix density: ", dens)
		if(dens >= 0){
            if( dens == 0 ){
				rhofac = 0.5
            }else{
                if( rhofac != 0 & prev_dens > 0 && dens>0 && dens != prev_dens ){
					p=log( iDens / dens ) / log( dens / prev_dens )
                    rhofac = rhofac^p
					if(rhofac>10){
						rhofac=10
					}
                }else{
					if(iDens > dens){
						rhofac = 0.9
					}else{rhofac = 1.1}
                }
            }
			rho = rho*rhofac
		}
		if(rho<=0 | rho>=1 | I==nIter){
			message("estimateIvcov() matrix density: ", bestDense)
			message("estimateIvcov() rho: ", rho)
			message("estimateIvcov() BestRho: ", bestRho)
			rho=bestRho
			rhoMtx=mkRhoMatrix(d,dis_idx,rho,BIG)
			a=glasso(X,rho=rhoMtx,zero=dis_idx)
#a=glasso(X,rho=rho)
			wi=a$wi
			return(wi)
		}
		message("estimateIvcov() rho: ", rho)
		rhoMtx=mkRhoMatrix(d,dis_idx,rho,BIG)
		a=glasso(X,rho=rhoMtx,zero=dis_idx)
#a=glasso(X,rho=rho)
		wi=a$wi
		prev_dens=dens
		dens=matrixDensity(wi)
		denseDiff=abs(iDens - dens)/dens
		if(denseDiff<bestDenseDiff){
			bestDenseDiff=denseDiff
			bestRho=rho
			bestDense=dens
		}
		I=I+1
	}
	wi
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
	}else if(opt$transform>1){
		write("Error: Unknown input matrix transformation mode",stderr())
		q(status=1)
	}
	X=estimateIvcov(mkSymmetric(C),dens)
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

