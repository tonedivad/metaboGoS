### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Return MS/MS sparse matrix from a list spectra
#' 
#' @param speclist List of spectra matrices
#' @param bin m/z binning
#' @param binwin  m/z window
#' @param reso  MS/MS resolution 
#' @param intlim Limit on the intensity normalised to the most intense peak
#' @param intlim2 Limit on the intensity as raw
#' @param nmax number of rows in the sparse matrix, i.e. max(mz)/bin
#' @import Matrix
#' 
#' @export
combMSMS<-function(speclist,bin=0.01,binwin=0.1,reso=17500,intlim=10^-5,intlim2=200,nmax=NA,mzlim=c(-Inf,Inf)){
  
  # bin=0.01;binwin=0.1;reso=17500;intlim=10^-5;intlim2=200;nmax=NA;mzlim=c(50,2000)
  
  avect=vector(mode = "list",length=length(speclist))
  names(avect)=names(speclist)
  for(im in 1:length(speclist)){
    if(im%%100==0) cat(ifelse(im%%1000==0,"X",'x'))
    m=speclist[[im]]
    sp2=do.call("rbind",lapply(1:nrow(m),function(i){
      m0=m[i,"mz"]
      gamma =m0/reso
      v=seq(m0-binwin/2,m0+binwin/2,bin/10)
      v=tapply((gamma^2)/(gamma^2 + (v-m0)^2),round(v/bin),sum)
      v=v*m[i,"y"]/sum(v)
      cbind(as.numeric(names(v)),v)
    }))
    sp2=tapply(sp2[,2],sp2[,1],sum)
    sp2=cbind(mz=as.numeric(names(sp2)),y=sp2)
    sp2=sp2[which((sp2[,"y"]/max(sp2[,"y"]))>=intlim & sp2[,"y"]>=intlim2),,drop=F]
    sp2=sp2[which((sp2[,"mz"]*bin)>=mzlim[1] & (sp2[,"mz"]*bin)<=mzlim[2]),,drop=F]
    if(nrow(sp2)==0) next
    avect[[im]]=sp2
  }

  if(is.na(nmax)) nmax=ceiling(max(sapply(avect,function(x) range(x[,"mz"]))))
  
  
  msp=do.call(.MGsvcbind, lapply(avect,function(x) sparseVector(i = x[,"mz"],x=x[,"y"],length = nmax)))
  colnames(msp)=names(avect)
  spform=paste0("%.",ceiling(abs(log10(bin)))+1,"f")
  rownames(msp)=sprintf(spform,(1:nmax)*bin)
  msp=msp[which(rowSums(msp>0)>0)[1]:rev(which(rowSums(msp>0)>0))[1],,drop=F]
  invisible(msp)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Compute smilarity between 1 or 2 MS/MS sparse matrices
#' 
#' @param m1 Sparse matrix
#' @param m2 Sparse
#' @param weight  Weight for NDP c(m=0.5,n=1.5)
#' @param plim set anything below plim to zero
#' @param useSqrt should sqrt tranform the intensities
#' @import Matrix
#' 
#' @export
compSimSPmat<-function(m1,m2=NULL,weight=c(m=NA,n=NA),plim=10^-5,useSqrt=TRUE){
  
  # weight=c(m=NA,n=NA);plim=10^-5;coslim=0;useSqrt=TRUE
  
  if(!is.null(m2)){
    
    same=(nrow(m1)==nrow(m2))
    if(same) same=(all(rownames(m1)==rownames(m2)))
    
    if(!same){
      cat("Two matrices are of different sizes\n")
      l1=names(which(rowSums(m1!=0)>0))
      l2=names(which(rowSums(m2!=0)>0))
      lint=intersect(l1,l2)
      l22=l1[!l1%in%lint]
      l21=l2[!l2%in%lint]
      m1=rbind(m1[l1,,drop=FALSE],Matrix(matrix(0,nrow=length(l21),ncol=ncol(m1),dimnames = list(l21,colnames(m1))), sparse = TRUE))
      m2=rbind(m2[l2,,drop=FALSE],Matrix(matrix(0,nrow=length(l22),ncol=ncol(m2),dimnames = list(l22,colnames(m2))), sparse = TRUE))
      lso=lso=order(as.numeric(rownames(m1)))
      m1=m1[lso,,drop=FALSE]
      m2=m2[rownames(m1),,drop=FALSE]
      
    }
  } else m2=m1=m1[which(rowSums(m1!=0)>0),,drop=FALSE]
  
  ## normalise
  m1=sweep(m1,2,apply(m1,2,max),"/")
  m1[which(m1<plim)]=0
  m2=sweep(m2,2,apply(m2,2,max),"/")
  m2[which(m2<plim)]=0
 
  if(useSqrt){ 
  m1=sqrt(m1)
  m2=sqrt(m2)
}
  ##
  percexp=sweep(crossprod(m1,m2>0),1,colSums(m1),"/")*100
  ## Comp NDP
  if(!all(is.na(weight))){
    mz=as.numeric(rownames(m1))^weight["m"]
  m1=sweep(m1^weight["n"],1,mz,"*")
  m2=sweep(m2^weight["n"],1,mz,"*")
  }
  
  ## Comp similarity
  cosine=crossprod(m1,m2)/sqrt(colSums(m1^2))%*%t(sqrt(colSums(m2^2)))

  invisible(list(cosine=cosine,percexp=percexp))
  
  
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Faster cbind function(for list of sparse vector)
#' see http://stackoverflow.com/questions/8843700/creating-sparse-matrix-from-a-list-of-sparse-vectors
#' 
#' @param ... whatever
#' @import Matrix
#' @return sparse matrix
#' @keywords internal
#' 
#' @export
.MGsvcbind <- function (...) {
  input <- lapply( list(...), as, "dsparseVector" )
  thelength <- unique(sapply(input,length))
  stopifnot( length(thelength)==1 )
  return( sparseMatrix( 
    x=unlist(lapply(input,slot,"x")), 
    i=unlist(lapply(input,slot,"i")), 
    p=c(0,cumsum(sapply(input,function(x){length(x@x)}))),
    dims=c(thelength,length(input))
  ) )
}

