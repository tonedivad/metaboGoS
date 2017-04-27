
# compSimSpectra<-function(thspec,exspec,mztol=0.1,nMatch=3,lim=0.01,xrankScoremat,verbose=T){
#   
#   .prepspec<-function(aspec,lim=0.01,typ="L"){
#     id=1:nrow(aspec)
#     if(any(c("x","mz")%in%colnames(aspec))) imz=aspec[,which(colnames(aspec)%in%c("x","mz"))[1]] else imz=aspec[,1]
#     if(any(c("y","ynorm","int","intnorm")%in%colnames(aspec))) iint=aspec[,which(colnames(aspec)%in%c("y","ynorm","int","intnorm"))[1]] else iint=aspec[,2]
#     m=cbind(mz=imz,int=iint,id=id,rk=1)
#     m[,2]=round(100*m[,2]/max(m[,2]),ceiling(-log10(0.01))+1)
#     m=m[order(m[,1]),,drop=F]
#     m=m[which(m[,2]>lim),,drop=F]
#     m[,"rk"]=rank(-m[,2],ties="max")
#     colnames(m)=paste0(colnames(m),typ)
#     m
#   }
#   
#   nthspec=lapply(thspec,.prepspec,lim,"L")
#   if(is.null(names(nthspec))) names(nthspec)=paste0("Lib",1:nrow(nthspec))
#   nexspec=lapply(exspec,.prepspec,lim,"U")
#   if(is.null(names(nexspec))) names(nexspec)=paste0("Exp",1:nrow(nexspec))
#   
#   
#   aresMatch=NULL
#   ares=list()
#   for(iexp in names(nexspec)){
#     if(verbose) cat(" ** ",iexp," ",sep="")
#     for(x in names(thspec)){
#       if(verbose) cat(".")
#       re=.GR2compSimCosM(munk = nexspec[[iexp]],mlib = nthspec[[x]],mztol = mztol,nMatch = nMatch,xrankScoremat=xrankScoremat,idunk=iexp,idlib=x)
#       if(is.null(re)) if(verbose) cat(" pb with ",x," ",sep="")
#       if(!is.null(re)) ares=c(ares,list(re))
#     }
#     if(verbose) cat("\n",sep="")
#   }
#   ares=ares[which(sapply(ares,length)==2)]
#   if(length(ares)>0){
#     aresMatch=do.call("rbind",lapply(ares,function(x) x$Stats))
#     rownames(aresMatch)=NULL
#     aresMatch=aresMatch[order(-aresMatch$Cosine),]
#   }
#   return(aresMatch)
#   
# }


# .GR2matchpairs<-function(munk,mlib,mztol=0.1,nMatch=3){
#   
#   ##### Get the best one to one match between spectra
#   system.time(resma<-.MGmatchtwosets(munk[,1],mlib[,1],tol = mztol))
#   if(is.null(resma))  return(NULL)
#   
#   resma=cbind(resma,munk[resma[,1],],mlib[resma[,2],])
#   resma=cbind(resma,scU=2*log(resma[,4])+.5*log(resma[,5]),
#               scL=2*log(resma[,6]+.5*log(resma[,7])))
#   ## ties brocken according to MassBank paper
#   resma=cbind(resma,Sc=resma[,"scU"]+resma[,"scL"])
#   if(any(resma[,3]>0)){
#     l=which(resma[,3]>0)
#     system.time(l2k<-.GR2inbip(e1=paste('U',resma[l,1],sep=''),e2=paste('L',resma[l,2],sep=''),we=resma[l,'Sc'])) 
#     if(length(l2k)!=length(l)) resma=resma[-l[-l2k],,drop=FALSE]
#   }
#   l2rm=c(which(duplicated(resma[,1]) & !is.na(resma[,1])),which(duplicated(resma[,2]) & !is.na(resma[,2])))
#   if(length(l2rm)>0) resma=resma[-unique(l2rm),,drop=F]
# 
#   resma
# }

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' cosine similarity
#'
#' @param munk unknown spectra 
#' @param mlib library spectra 
#' @param mztol tolerance
#' @param nMatch min number of matches
#' @param xrankScoremat xrankScoremat
#' @param idunk id of unknown
#' @param idlib id of library spectra
#' @return matrix of matches
#' @export
compOneSimCosM<-function(munk,mlib,mztol=0.1,nMatch=3,xrankScoremat=NULL,idunk="unk",idlib="lib"){
  
  ##### Get the best one to one match between spectra
  system.time(resma<-.MGmatchtwosets(munk[,1],mlib[,1],tol = mztol))
  if(is.null(resma))  return(NULL)
  colnames(resma)=c("row","col","typ")
  
  munk[,2]=100*munk[,2]/max(munk[,2])
  colnames(munk)[1:2]=c("mzU","intU")
  mlib[,2]=100*mlib[,2]/max(mlib[,2])
  colnames(mlib)[1:2]=c("mzL","intL")
  if(!is.null(xrankScoremat)){
    munk=cbind(munk,"rkU"=rank(-m[,'intU'],ties="max"))
    mlib=cbind(mlib,"rkL"=rank(-m[,'intL'],ties="max"))
  }
  
  
  resma=cbind(resma,munk[resma[,1],],mlib[resma[,2],])
  resma=cbind(resma,scU=2*log(resma[,'mzU'])+.5*log(resma[,'intU']),
              scL=2*log(resma[,'mzL']+.5*log(resma[,'intL'])))
  ## ties brocken according to MassBank paper
  resma=cbind(resma,Sc=resma[,"scU"]+resma[,"scL"])
  if(any(resma[,3]>0)){
    l=which(resma[,3]>0)
    system.time(l2k<-.GR2inbip(e1=paste('U',resma[l,1],sep=''),e2=paste('L',resma[l,2],sep=''),we=resma[l,'Sc'])) 
    if(length(l2k)!=length(l)) resma=resma[-l[-l2k],,drop=FALSE]
  }
  l2rm=c(which(duplicated(resma[,1]) & !is.na(resma[,1])),which(duplicated(resma[,2]) & !is.na(resma[,2])))
  if(length(l2rm)>0) resma=resma[-unique(l2rm),,drop=F]
  
  resma=resma[order(apply(resma[,c('mzU','mzL')],1,median,na.rm=T)),,drop=F]
  
  ##### 
  Stats=c(Nu=sum(!is.na(resma[,"mzU"])),Nl=sum(!is.na(resma[,"mzL"])),Ncom=sum(!is.na(resma[,"mzU"]) & !is.na(resma[,"mzL"])),
          Cosine=NA,XScore=NA,XScoreU=NA,Pearson=NA)
  
  
  ## Stats
  if(Stats["Ncom"]>=nMatch){
    Stats["Cosine"]=round(sum((resma[,"intU"]*resma[,"intL"])/sqrt(sum(resma[,"intU"]^2,na.rm=T)*sum(resma[,"intL"]^2,na.rm=T)),na.rm=T),4)
    Stats["Pearson"]=round(cor.test(resma[,"intU"],resma[,"intL"],method="pearson")$estimate,4)
    if(!is.null(xrankScoremat)){
      lrkU=which(resma[,"rkU"]<31)
      mU=resma[lrkU,c("rkL","rkU")]
      mU[which(mU[,1]>30 | is.na(mU[,1])),1]=31
      xmu=mean(xrankScoremat[mU])
      lrkL=which(resma[,"rkL"]<31)
      mL=resma[lrkL,c("rkU","rkL")]
      mL[which(mL[,1]>30 | is.na(mL[,1])),1]=31
      xml=mean(xrankScoremat[mL])
      Stats["XScoreU"]=xmu
      Stats["XScore"]=(xml+xmu)/2
    }
    
  }
  Stats=data.frame(Unk=idunk,Lib=idlib,t(Stats))
  res=list(Match=resma,Stats=Stats)
  class(res)<-c(class(res),'compOneSimCosM')
  return(res)
}

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' find best pairs
#'
#' @param v1 vector of m/z 
#' @param v2 vector of m/z 
#' @param tol tolerance
#' @param isdppm tolerance in ppm
#' @param long long matrix
#' @keywords internal
#' @return matrix of matches
#' @export

.MGmatchtwosets<-function(v1,v2,tol=1,isdppm=FALSE,long=TRUE){
  if(!isdppm) mma=which(abs(outer(v1,v2, "-"))<=tol,arr.ind = T) else mma=which(abs(outer(v1,v2, "/")-1)*10^6<=tol,arr.ind = T)
  if(nrow(mma)==0) return(NULL)
  mma=cbind(mma,mma[,1]%in%unique(mma[duplicated(mma[,1]),1]),mma[,2]%in%unique(mma[duplicated(mma[,2]),2]))
  mma1=rbind(cbind(1:length(v1),NA,-1),cbind(NA,1:length(v2),-1))
  # mma1=mma1[which(!(mma1[,1]%in%mma[,1] |  mma1[,2]%in%mma[,2])),,drop=F]
  if(long) return(rbind(cbind(mma[,1:2,drop=FALSE],rowSums(mma[,3:4,drop=F])),mma1))
  #### Case: no dups
  mma0=NULL
  l=which(rowSums(mma[,3:4,drop=F])==0)
  if(length(l)>0){
    mma0=rbind(mma0,cbind(mma[l,1:2,drop=FALSE],0,(v1[mma[l,1]]-v2[mma[l,2]])))
    mma=mma[-l,,drop=F]
  }
  if(nrow(mma)==0) return(mma0)
  #### Case: dups in 1->unique in 2
  l=which(mma[,1]%in%names(which(tapply(mma[,4],mma[,1],sum)==0)))
  if(length(l)>0){
    mma0=rbind(mma0,cbind(mma[l,1:2,drop=F],1,(v1[mma[l,1]]-v2[mma[l,2]])))
    mma=mma[-l,,drop=F]
  }
  if(nrow(mma)==0) return(mma0)
  #### Case: dups in 2->unique in 1
  l=which(mma[,2]%in%names(which(tapply(mma[,3],mma[,2],sum)==0)))
  if(length(l)>0){
    mma0=rbind(mma0,cbind(mma[l,1:2,drop=F],2,(v1[mma[l,1]]-v2[mma[l,2]])))
    mma=mma[-l,,drop=F]
  }
  if(nrow(mma)==0) return(mma0)
  mma0=rbind(mma0,cbind(mma[,1:2,drop=FALSE],3,(v1[mma[,1]]-v2[mma[,2]])))
  return(mma0)
}

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' weighted bipartite matching
#'
#' @param e1 list of edges 
#' @param e2 list of edges 
#' @param we weight
#' @import igraph
#' @keywords internal
#' @return matrix of matches
#' @export

.GR2inbip<-function(e1,e2,we){
  library(igraph)
  ig=igraph:::graph_from_edgelist(cbind(e1,e2), directed = FALSE)
  V(ig)$type <- grepl('U',V(ig)$name)
  retmp=max_bipartite_match(ig,weights =we)$matching
  l=which(paste(e1,e2)%in%paste(retmp,names(retmp)))
  return(l)
}

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' plotting 
#'
#' @param resMatch result from compOneSimCosM match 
#' @param typ 1/2
#' @import maptools
#' @return matrix of matches
#' @export

plotResMatch<-function(resMatch,type=1){
  
  if(!"compOneSimCosM"%in%class(resMatch)) stop('Not the right class')
  resma=resMatch[[1]]
  yaxt=seq(-100,100,20)
  xaxt=pretty(c(resma[,"mzU"],resma[,"mzL"]))
  
  ## stats
  perc=100*colSums(resma[which(resma[,"typ"]==0),c("intU","intL")])/colSums(resma[,c("intU","intL")],na.rm=T)
  n=colSums(!is.na(resma[,c("intU","intL")]))
  main=sprintf("Common:%d, %s(%d, %.2f), %s(%d, %.2f)",resMatch$Stats$Ncom,resMatch$Stats$Unk,n[1],perc[1],resMatch$Stats$Lib,n[2],perc[2])
  smain=sprintf("Cosine=%.3f, PCC=%.2f, xScore=%.3f",resMatch$Stats$Cosine,resMatch$Stats$Pearson,resMatch$Stats$XScore)
  
  
  if(type==1){
    plot(range(xaxt),range(yaxt),xlim = range(xaxt),ylim = range(yaxt)*1.1,cex=0,axes=F,xlab="",ylab="Rel int.")
    title(smain,sub=main,cex.main=1)
    axis(1,at=xaxt,pos=0,cex.axis=.5)
    axis(2,at=yaxt,las=2)
    l=which(is.na(resma[,1]));if(length(l)>0) segments(resma[l,"mzL"],0,resma[l,"mzL"],-resma[l,"intL"],col="grey20")
    l=which(is.na(resma[,2]));if(length(l)>0) segments(resma[l,"mzU"],0,resma[l,"mzU"],resma[l,"intU"],col="grey20")
    l=which(!is.na(resma[,1]) & !is.na(resma[,2]))
    if(length(l)>0) segments(resma[l,"mzL"],0,resma[l,"mzL"],-resma[l,"intL"],col="red",lwd=2)
    if(length(l)>0) segments(resma[l,"mzU"],0,resma[l,"mzU"],resma[l,"intU"],col="blue",lwd=2)
    l=which(!is.na(resma[,1]) & !is.na(resma[,2]) & resma[,"intU"]>5)
    if(length(l)>0){
      labs=round(resma[l,"mzU"],3)
      pointLabel(resma[l,"mzU"],resma[l,"intU"],labels = as.character(labs))
    }
    l=which(!is.na(resma[,1]) & !is.na(resma[,2]) & resma[,"intL"]>5)
    if(length(l)>0){
      labs=round(resma[l,"mzL"],3)
      pointLabel(resma[l,"mzL"],-resma[l,"intL"],labels = as.character(labs))
    }
  }
  
  if(type==2){
    yaxt=seq(0,100,10)
    l=which(!is.na(resma[,1]) & !is.na(resma[,2]))
    plot(resma[l,"intL"],resma[l,"intU"],axes=F,cex=1,pch=16,xlim = range(yaxt)*1.1,ylim = range(yaxt)*1.1,xlab="Library int.",ylab="Unknown int.")
    title(smain,sub=main,cex.main=1)
    segments(-3,-3,103,103)
    axis(1,at=yaxt,pos=-3);axis(2,at=yaxt,las=2,pos=-3)
    labs=round(apply(resma[l,c('mzU','mzL')],1,median,na.rm=T),3)
    pointLabel(resma[l,"intL"],resma[l,"intU"],labels = as.character(labs))
  }
}


