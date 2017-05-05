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
#' @param reduce should all zeros row removed
#' @import Matrix
#' 
#' @export
combMSMS<-function(speclist,bin=0.01,binwin=0.1,reso=17500,intlim=10^-5,intlim2=200,nmax=NA,mzlim=c(-Inf,Inf),reduce=F){
  
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
  msp=msp[which(Matrix:::rowSums(msp)>0)[1]:rev(which(Matrix:::rowSums(msp)>0))[1],,drop=F]
  if(reduce) msp=msp[which(Matrix:::rowSums(msp)>0),,drop=F]
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
#' @param retPerc return perc overlap
#' @import Matrix
#' 
#' @export
compSimSPmat<-function(m1,m2=NULL,weight=c(m=NA,n=NA),plim=10^-5,useSqrt=TRUE,retPerc=TRUE){
  
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
  } else m2=m1=m1[which(Matrix:::rowSums(m1!=0)>0),,drop=FALSE]
  
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
  if(retPerc) percexp=sweep(Matrix:::crossprod(m1,m2>0),1,Matrix:::colSums(m1),"/")*100
  ## Comp NDP
  if(!all(is.na(weight))){
    mz=as.numeric(rownames(m1))^weight["m"]
    m1=sweep(m1^weight["n"],1,mz,"*")
    m2=sweep(m2^weight["n"],1,mz,"*")
  }
  
  ## Comp similarity
  cosine=crossprod(m1,m2)/sqrt(colSums(m1^2))%*%t(sqrt(colSums(m2^2)))
  
  if(retPerc) return(list(cosine=cosine,percexp=percexp)) else return(cosine)
  
  
  
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



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Merge list of ms/ms spectra
#' 
#' @param split list of spectra
#' @param dmz delta m/z for merging
#' @param fct  merging function
#' @param minNoise supress anything below minNoise
#' 
#' @export
mergeMSMS<-function(splist,dmz=0.005,fct=median,minNoise=0,reso=17500){
  
  whichmz="mz"
  whichint="y"
  whichid="sc"
  xtmp=do.call("rbind",splist)[,c(whichid,whichmz,whichint)]
  if(nrow(xtmp)<2) return(xtmp)
  xtmp=xtmp[order(xtmp[,whichmz],xtmp[,whichint]),,drop=F]
  llsplit=lapply(.GRsplistMZ(xtmp[,"mz"],dmz=dmz*4.1),function(x) xtmp[x,,drop=F])
  
  for(i in which(sapply(llsplit,nrow)>1)){
    itmp=llsplit[[i]]
    if(length(unique(itmp[,1]))==1) next
    ll=tapply(1:nrow(itmp),itmp[,whichid],function(y) itmp[y,c(whichmz,whichint),drop=F])
    cspec=combMSMS(ll,bin = dmz/10,reso = reso,reduce = F,mzlim = range(itmp[,whichmz])+c(-2,2)*dmz,intlim = 0,intlim2 = 0)
    cspec=cbind(as.numeric(rownames(cspec)),rowSums(cspec))
    nsp=2
    lpks=.GRfindturnpoint(cspec[,2])$pks
    if(length(lpks)==0){nsp=1;lpks=which(.GRipeaks(-diff(diff(cspec[,2])),span = 2*nsp+1))+1}
    lpks=outer(lpks,-nsp:nsp,"+")           
    lpks[which(lpks<1 | lpks>nrow(cspec))]=NA
    newsp=t(apply(lpks,1,function(x) c(mz=weighted.mean(cspec[x,1],cspec[x,2],na.rm = T),y=max(cspec[x,2],na.rm = T))))
    newsp[,2]=newsp[,2]*sum(itmp[,whichint])/sum(newsp[,2])
    llsplit[[i]]=cbind(0,newsp)
  }
  finalcs=do.call("rbind",llsplit)
  colnames(finalcs)=c(whichid,whichmz,whichint)
  finalcs[,3]= finalcs[,3]*sum(xtmp[,whichint])/sum( finalcs[,3])
  finalcs[order(finalcs[,2]),,drop=F]
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Merge list of ms/ms spectra
#' 
#' @param split list of spectra
#' @param dmz delta m/z for merging
#' @param fct  merging function
#' 
#' @export
mergeMSMS.2<-function(splist,dmz=0.001,fct=sum){
  whichmz="mz"
  whichint="y"
  nspec=length(splist)
  xtmp=do.call("rbind",lapply(1:nspec,function(x) cbind(sc=x,mz=splist[[x]][,whichmz],y=splist[[x]][,whichint])))
  if(nrow(xtmp)<2) return(xtmp)
  xtmp=xtmp[order(xtmp[,whichmz],xtmp[,whichint]),,drop=F]
  xtmp=xtmp[order(xtmp[,whichmz],xtmp[,whichint]),,drop=F]
  diffmz=diff(xtmp[,whichmz])
  while(any(diffmz<dmz) & nrow(xtmp)>1){
    (imerg=which.min(diffmz))
    # if(!keepMax){
    xtmp[imerg+1,whichmz]=weighted.mean(xtmp[imerg+(0:1),whichmz],xtmp[imerg+(0:1),whichint])
    xtmp[imerg+1,whichint]=sum(xtmp[imerg+(0:1),whichint])
    #  } else imerg=imerg+which.min(xtmp[imerg+(0:1),whichint])-1
    xtmp=xtmp[-imerg,]
    xtmp=xtmp[order(xtmp[,whichmz]),,drop=F]
    diffmz=diff(xtmp[,whichmz])
  }
  xtmp
}