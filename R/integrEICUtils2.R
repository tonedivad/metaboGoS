### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name .MGmultiEICintgrPar2
#' @title multiple integration of several eic
#' 
#' Multicore wrapper for .MGmultiEICintgr
#'
#' @param alleicmat eic: list of mz, y, and y2
#' @param parDeco deconvolution parameter list
#' @param l2excl vector of sample ids to be excluded when refining peak allocations
#' @param doPlot plotting - disabled if run in parallele
#' @param nSlaves num slaves
#' @param minPeakHeight minPeakHeight
#' @param letSa sample id
#' @param colSa sample colors
#' @param typSa sample line types
#' @import matrixStats
#' @return peak matrix data.frame
#' @export
.MGmultiEICintgrPar2<-function(alleicmat,parDeco,l2excl=NULL,doPlot=TRUE,nSlaves=1,minPeakHeight=NULL,letSa=NULL,colSa=NULL,typSa=NULL){
  
  
  # l2excl=NULL;doPlot=1;nSlaves=1;minPeakHeight=NULL
  
  if(is.null(minPeakHeight)) minPeakHeight=max(parDeco$minHeightMS1,parDeco$minNoiseMS1*parDeco$sbr)
  
  llre=list()
  ll=names(alleicmat)#[111:120]
  ll=ll[which(sapply(alleicmat[ll],function(x) max(x$y,na.rm=T))>=minPeakHeight)]
  lperc=""
  if(length(ll)>20) lperc=ll[round(seq(1,length(ll),length=12)[2:11])]
  
  
  if(nSlaves>1)   nSlaves=max(1, min(nSlaves,detectCores()-1))
  if(nSlaves>1){
    clProc<-makeCluster(nSlaves)
    doParallel::registerDoParallel(clProc)
    cat(" -- registering ",nSlaves," clusters\n",sep="")
  }
  
  
  ## Fix color/letters
  if(length(doPlot) & nSlaves==1){
    lsids=sort(unique(unlist(sapply(alleicmat,function(x) rownames(x[[1]])))))
    if(is.null(letSa)){
      letSa=rep(c(letters,LETTERS),ceiling(length(lsids)/52))[1:length(lsids)]
      names(letSa)=lsids
    }
    if(is.null(colSa)){
      colSa=rep(brewer.pal(8,"Dark2"),ceiling(length(lsids)/8))[1:length(lsids)]
      names(colSa)=lsids
    }
    if(is.null(typSa)){
      typSa=(lsids%in%l2excl)+1
      names(typSa)=lsids
    }
  }
  
  ### Parallele bit
  if(nSlaves>1)
    llre=foreach(idx = ll,.packages = c("metaboGoS"), .verbose =FALSE)  %dopar%{
      re=.MGmultiEICintgr2(eic=alleicmat[[idx]],parDeco,doPlot = FALSE,l2excl = l2excl,minPeakHeight=minPeakHeight)
      if(is.null(re)) return(idx)
      list(data.frame(RoiId=idx,re$MatPks),data.frame(RoiId=idx,re$Pk))
    }
  ## Serial bit
  if(nSlaves<=1) for(idx in ll){
    # if(nSlaves<=1) for(idx in sample(ll[!ll%in%names(llre)])){
    if(idx %in% lperc) cat(idx,"(",which(ll==idx),") ",sep="")
    re=.MGmultiEICintgr2(eic=alleicmat[[idx]],parDeco,doPlot = doPlot,l2excl = l2excl,minPeakHeight=minPeakHeight,colSa = colSa,letSa=letSa,typSa = typSa,main = idx)
    if(is.null(re)) llre[[idx]]=idx else  llre[[idx]]=list(data.frame(RoiId=idx,re$MatPks),data.frame(RoiId=idx,re$Pk))
  }
  
  if(nSlaves>1) stopCluster(clProc)
  
  ## combine results
  lndf=which(!sapply(llre,is.list))
  if(length(lndf)) cat(" ++++ ",length(lndf)," excluded ROIs:\n",paste(names(lndf),collapse = " "),"\n",sep="")
  llre=llre[which(sapply(llre,is.list))]
  if(length(llre)==0){
    cat("Something went wrong no acceptable ROI were found!!!")
    return(NULL)
  }
  amatpks=do.call("rbind",lapply(llre,function(x) x[[1]]))
  apks=do.call("rbind",lapply(llre,function(x) x[[2]]))
  
  return(list(MatPks=amatpks,Pks=apks))
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name .MGmultiEICintgr2
#' @title multiple integration of several eic
#'
#' @param eic eic: list of mz, y, and y2
#' @param parDeco deconvolution parameter list
#' @param l2excl vector of sample ids to be excluded when refining peak allocations
#' @param doPlot plotting at each step, vector of 0,1,2
#' @param minPeakHeight minPeakHeight
#' @param letSa sample id
#' @param colSa sample colors
#' @param typSa sample line types
#' @param main title
#' @import matrixStats
#' @import NMF
#' @return peak matrix data.frame
#' @export
.MGmultiEICintgr2<-function(eic,parDeco=list(psdrt=0.00446,span=7,
                                             minHeightMS1=100000,minNoiseMS1=5000),l2excl=NULL,doPlot=0,minPeakHeight=10^6,
                            colSa=NULL,letSa=NULL,typSa=NULL,main=NULL){
  

  ########## ############ ########## ############ ##########
  ## Prep data 
  ########## ############ ########## ############ ##########
  m=t(eic$y2)#/1000000
  m0=t(eic$y)#/1000000
  mmz=t(eic$mz)#/1000000
  if(all(colnames(m)%in%l2excl)) return(NULL) ## if all to be excluded then stop
  maxm0=max(max(m0[,!colnames(m0)%in%l2excl],na.rm=T),max(m[,!colnames(m)%in%l2excl],na.rm=T))
  if(maxm0<minPeakHeight) return(NULL)
  
  span=parDeco$span
  minNoise=parDeco$minNoiseMS1
  nspan=floor(span/2)*2+1
  span2=ifelse(is.null(parDeco$Rois$span2),2*nspan+1,parDeco$Rois$span2)
  
  ## Fix color/letters
  if(is.null(letSa) & length(doPlot)){
    letSa=rep(c(letters,LETTERS),ceiling(ncol(m)/52))[1:ncol(m)]
    names(letSa)=colnames(m)
  }
  if(is.null(colSa) & length(doPlot)){
    colSa=rep(brewer.pal(8,"Dark2"),ceiling(ncol(m)/8))[1:ncol(m)]
    names(colSa)=colnames(m)
  }
  if(is.null(typSa) & length(doPlot)){
    typSa=(colnames(m0)%in%l2excl)+1
    names(typSa)=colnames(m)
  }
  
  ########## ############ ########## ############ ##########
  ## Generate good peaks 
  ########## ############ ########## ############ ##########
  if(!all(c( "bsl","bslc")%in%names(eic))){
    mbsl=mbslsc=m
    for(isid in colnames(m)){
      re=.MGinfctonebsl(m[,isid],span=span,minNoise=minNoise)
      mbsl[,isid] = re[[1]]
      mbslsc[,isid] = re[[2]]
    }
  }else{
    mbsl=t(eic$bsl)
    mbslsc=t(eic$bslc)
  }
  ########
  apks=list()
  for(isid in colnames(m)){
    # x=1:nrow(m);y=m[,isid];bsl = mbsl[,isid];bslscore = mbslsc[,isid];snr.thresh =parDeco$sbsl;  span=nspan;minNoise = minNoise*1.01;v2alim = 0.8
    re=.MGsimpleIntegr2(x=1:nrow(m),y=m[,isid],bsl = mbsl[,isid],bslscore = mbslsc[,isid],snr.thresh =parDeco$sbslr,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8,span2 = span2+2)
    if(nrow(re)==0) next
    re$Sid=isid
    apks[[isid]]=re
  }
  apks=do.call("rbind",apks)
  if(is.null(apks)){
    if(doPlot) .infctintegrplot(m,m0,parDeco,NULL,fac=10^6,rmz=range(mmz,na.rm=T),cols =colSa,sidsl = letSa,typs=typSa,main=main)
    return(NULL)
  }
  names(apks)=gsub("pk\\.","tick.",names(apks))
  apks$PkCl=apks$tick.max=NA
  ## get the max m0 in the peak region
  apks$tick.max=apply(apks,1,function(x) max(m0[x[2]:x[3],x['Sid']],na.rm=T))
  apks$Above=(apks$tick.int>=minPeakHeight | apks$tick.max>=minPeakHeight)
  if(all(!apks$Above)){
    if(doPlot) .infctintegrplot(m,m0,parDeco,NULL,fac=10^6,rmz=range(mmz,na.rm=T),cols =colSa,sidsl = letSa,typs=typSa,main=main)
    return(NULL)
  }
  
  
  ########## ############ ########## ############ ##########
  ## Compute master peak 
  ########## ############ ########## ############ ##########
  lcol2use=which(!colnames(m)%in%l2excl & colnames(m)%in%apks$Sid) ## before
  if(length(lcol2use)==0){
    if(doPlot) .infctintegrplot(m,m0,parDeco,NULL,fac=10^6,rmz=range(mmz,na.rm=T),cols =colSa,sidsl = letSa,typs=typSa,main=main)
    return(NULL)
  }
  if(length(lcol2use)>1){
    
    sv <- svd(m[,lcol2use,drop=F], nu = 1, nv = 1)
    renmfv=round(m[,lcol2use,drop=F]%*%sv$v,3)
    renmfvnoise=round(sum(parDeco$minNoiseMS1*sv$v),3)
    
  if(quantile(renmfv,.8)<0){
    renmfv=-renmfv
    renmfvnoise=-renmfvnoise
  }
  renmfv[renmfv<renmfvnoise]=renmfvnoise
  } else {renmfv=m[,lcol2use];renmfvnoise=parDeco$minNoiseMS1}
  ##### Compute baseline/peak bounds
  rebsl=.MGinfctonebsl(renmfv,span=span,minNoise=renmfvnoise)
  
 #  x=1:length(renmfv);y=renmfv;bsl=rebsl[[1]];bslscore=rebsl[[2]];minNoise = renmfvnoise;snr.thresh =1;span=nspan;v2alim = 0.8;span2 =span2
  
  pks=.MGsimpleIntegr2(1:length(renmfv),renmfv,bsl=rebsl[[1]] ,bslscore=rebsl[[2]],minNoise = renmfvnoise,snr.thresh =1,span=nspan,v2alim = 0.8,span2 =span2+2)
  if(nrow(pks)==0){
    if(doPlot) .infctintegrplot(m,m0,parDeco,NULL,fac=10^6,rmz=range(mmz,na.rm=T),cols =colSa,sidsl = letSa,typs=typSa,main=main,v=renmfv)
    return(NULL)
  }
  
  # lcol2use=which(!colnames(m)%in%l2excl)
  # if(length(lcol2use)>1){
  #   nDim=ifelse(length(lcol2use)>2,2,1)
  #   renmf<-try(NMF:::nmf(m[,lcol2use,drop=F], nDim, 'snmf/l',seed='ica'),T)
  #   while("try-error"%in%class(renmf) & nDim>1){
  #     nDim=nDim-1
  #     renmf<-try(NMF:::nmf(m[,lcol2use,drop=F], nDim, 'snmf/l',seed='ica'),T)
  #   }
  #   if(!"try-error"%in%class(renmf)) renmfv=rowSums(NMF:::loadings(renmf))
  #   if("try-error"%in%class(renmf)){
  #     cat("Switch to svd \n",sep="")
  #     renmfv=prcomp(m[,lcol2use,drop=F],center = FALSE)$x[,1]
  #     if(cor(rowSums(m[,lcol2use,drop=F]),renmfv)<0) renmfv=-renmfv
  #     renmfv[renmfv<0]=0
  #   }
  # } else renmfv=m[,lcol2use] 
  # if(nrow(pks)==0){
  #   cat(" -> snr down to 1 in ",main,"\n",sep="")
  #   pks=.MGsimpleIntegr2(1:length(renmfv),renmfv,bsl=rebsl[[1]] ,bslscore=rebsl[[2]],snr.thresh =.99,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8,span2 = 2*nspan+1)
  #   
  # }
  

  ### Only keep with high intensities clusters
  l2k=which(tapply(1:nrow(pks),pks$pk.cl,function(x){
    l=min(pks[x,1:3]):max(pks[x,1:3])
    max(max(m0[,lcol2use],na.rm=T),max(m[,lcol2use],na.rm=T))
  })>=minPeakHeight)
  pks=pks[pks$pk.cl%in%names(l2k),]
  if(nrow(pks)==0){
    if(doPlot) .infctintegrplot(m,m0,parDeco,NULL,fac=10^6,rmz=range(mmz,na.rm=T),cols =colSa,sidsl = letSa,typs=typSa,main=main,v=renmfv)
    return(NULL)
  }
  
  ### relabel
  pks$Pk=1:nrow(pks)
  pks$PkCl=as.numeric(factor(pks$pk.cl))
  

  
  # minpk=matrix(FALSE,ncol=ncol(m),nrow=nrow(m),dimnames = dimnames(m))
  # for(i in which(!apks$Sid%in%l2excl)) minpk[apks$tick.left[i]:apks$tick.right[i],apks$Sid[i]]=T
  # lcol2use=which(colSums(minpk)>0)
  # m2=m;m[!minpk]=0
  # renmfv=prcomp(m[,lcol2use,drop=F],center = FALSE)$x[,1]
  # if(cor(rowSums(m[,lcol2use,drop=F]),renmfv)<0) renmfv=-renmfv
  # 
  
  
  ########## ############ ########## ############ ##########
  ## Associate peak to clusters 
  ########## ############ ########## ############ ##########
  linpks=which(outer(apks$tick.loc,pks$pk.left,"-")>0 & outer(apks$tick.loc,pks$pk.right,"-")<0,arr.ind = T)
  if(nrow(linpks)==0){
    if(doPlot) .infctintegrplot(m,m0,parDeco,NULL,fac=10^6,rmz=range(mmz,na.rm=T),cols =colSa,sidsl = letSa,typs=typSa,main=main,v=renmfv)
    return(NULL)
  }
  
  if(max(table(linpks[,1]))>1) cat("Dups found in ",main,"\n",sep="")
  apks$PkCl[linpks[,1]]=linpks[,2]
  
  ## apksrm: keep the 'good' peaks somewhere
  apksrm=apks[1,][-1,]
  ## remove peaks that are NA/not above in their cluster
  lcl=paste0(apks$Sid,";;",apks$tick.cl)
  lcl2rm=names(which(tapply(!apks$Above & is.na(apks$PkCl),lcl,all)))
  if(length(lcl2rm)){
    apksrm=apks[which(lcl%in%lcl2rm),]
    apks=apks[which(!lcl%in%lcl2rm),]
    lcl=lcl[which(!lcl%in%lcl2rm)]
  }
  
  ### check if peaks that are NA/not above in their cluster may be left/right
  ## add contrain on apex hieght
  lim=nspan
  lfill=NULL
  linpksleft=which(outer(apks$tick.right,pks$pk.left,"-")>lim & outer(apks$tick.right,pks$pk.loc,"-")<0,arr.ind = T)
  if(nrow(linpksleft)) linpksleft=linpksleft[which(is.na(apks$PkCl[linpksleft[,1]])),,drop=F]
  if(nrow(linpksleft)){
    l2k=c()
    for(i in 1:(nrow(linpksleft))){
      x=linpksleft[i,]
      if(!is.na(apks$PkCl[x[1]])) next
      if(!x[2]%in%na.omit(apks$PkCl[lcl==lcl[x[1]]])) next
      if(mbslsc[apks$tick.right[x[1]],apks$Sid[x[1]]]<0) next
      l2k=c(l2k,i)
    }
    if(length(l2k)) lfill=rbind(lfill,linpksleft[l2k,,drop=F])
  }
  
  linpksright=which(outer(apks$tick.left,pks$pk.loc,"-")>0 & outer(apks$tick.left,pks$pk.right,"-")<(-lim),arr.ind = T)
  if(nrow(linpksright)) linpksright=linpksright[which(is.na(apks$PkCl[linpksright[,1]])),,drop=F]
  if(nrow(linpksright)){
    l2k=c()
    for(i in 1:(nrow(linpksright))){
      x=linpksright[i,]
      if(!x[2]%in%na.omit(apks$PkCl[lcl==lcl[x[1]]])) next
      if(mbslsc[apks$tick.left[x[1]],apks$Sid[x[1]]]<0) next
      l2k=c(l2k,i)
    }
    if(length(l2k)) lfill=rbind(lfill,linpksright[l2k,,drop=F])
  }
  
  if(!is.null(lfill)){
    lfill=lfill[!duplicated(lfill[,1]),,drop=F]
    if(nrow(lfill)) apks$PkCl[lfill[,1]]=lfill[,2]
  }
  
  ########## ########## ########## ##########
  # novel peak in l2excl
  l2rm=which(is.na(apks$PkCl) & apks$Sid%in%l2excl)
  if(length(l2rm)){
    apksrm=rbind(apksrm,apks[l2rm,])
    apks=apks[-l2rm,]
  }
  
  ########## ########## ########## ##########
  # remove peaks no above \ sbr<lim
  l2rm=which(is.na(apks$PkCl) & !(apks$Above & apks$tick.snr>=parDeco$sbslr))
  if(length(l2rm)){
    apksrm=rbind(apksrm,apks[l2rm,])
    apks=apks[-l2rm,]
  }
  

  ########## ########## ########## ##########
  # may have forgotten some good peak???
  lnew=which(is.na(apks$PkCl))
  if(length(lnew)){
    if(doPlot) cat("New peak found: ",length(lnew)," in ",main,"\n",sep="")
    apksrm=rbind(apksrm,apks[lnew,])
    apks=apks[-lnew,]
  }
  
  ########## ############ ########## ############ ##########
  ## Peak integration 
  ########## ############ ########## ############ ##########
  reint=.MGinEICintrfct(apks,m,m0,mbsl,mmz,parDeco,nspan)
  minpk=reint[[1]]
  matpks=reint[[2]]
  matpks$PkCl=pks$PkCl[matpks$Pk]
  matpks$InDeco=1 ## for post processing analysis
  matpks=matpks[order(matpks$PkCl,matpks$Pk,matpks$Sid,matpks$rt),]
  
  l=names(which(tapply(matpks$int.sm>=minPeakHeight | matpks$int.ap>=minPeakHeight,matpks$PkCl,any)))
#  print(pks$PkCl)
  
  if(!all(matpks$PkCl%in%l)){
    matpks$PkCl=match(matpks$PkCl,l)
    pks$PkCl=match(pks$PkCl,l)
    matpks=matpks[!is.na(matpks$PkCl),]
    pks=pks[!is.na(pks$PkCl),]
    
    
    ## make pk are from 1 to n
    lupk=unique(matpks$Pk)
    matpks$Pk=match(matpks$Pk,lupk)
    pks$Pk=match(pks$Pk,lupk)
    pks=pks[!is.na(pks$Pk),]
    
  }
  
  if(nrow(matpks)==0){
    if(doPlot) .infctintegrplot(m,m0,parDeco,NULL,fac=10^6,rmz=range(mmz,na.rm=T),cols =colSa,sidsl = letSa,typs=typSa,main=main,v=renmfv)
    return(NULL)
  }

  ########## ############ ########## ############ ##########
  ## Fill in missing 
  ########## ############ ########## ############ ##########
  tabmiss=table(factor(matpks$Sid,colnames(m)),factor(matpks$Pk))
  lmiss=which(tabmiss==0,arr=T)
  lmiss=data.frame(Sid=colnames(m)[lmiss[,1]],currpk=lmiss[,2])
  
  allre=list()
  if(nrow(lmiss)) for(i in 1:nrow(lmiss)){
    isid=lmiss$Sid[i]
    currpk=lmiss$currpk[i]
    ipkx=pks[pks$Pk==currpk,]
    lx=ipkx$pk.left:ipkx$pk.right
    v=m[lx,isid]
    v0=m0[,isid]
    if(all(is.na(v0[lx]))) next
    iv=minpk[lx,isid]
    l=which(matpks$Pk>currpk & matpks$Sid==isid);if(length(l)) iv[lx>min(matpks$rt[l])]=TRUE
    l=which(matpks$Pk<currpk & matpks$Sid==isid);if(length(l)) iv[lx<max(matpks$rt[l])]=TRUE
    if(all(v[!iv]<(minNoise*2))) next
    ## lower down the snr
    re=.MGsimpleIntegr2(x=lx,y=v,bsl = mbsl[lx,isid],bslscore = mbslsc[lx,isid],snr.thresh =1.01,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8,span2)
    if(nrow(re)){
      re=re[which.min(abs(re$pk.loc-ipkx$pk.loc)),]
      allre[[i]]=cbind(re,Sid=isid,tick.max=max(v0[re[1,2]:re[1,3]],na.rm=TRUE),PkCl=currpk)
      next
    }
    
    re=.MGsimpleIntegr(x=lx,y=v,noise.local = rep(minNoise,length(lx)) ,snr.thresh = 1.01,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8)
    if(nrow(re)){
      re=.MGinfctmergepk(re)
      re=re[which.min(abs(re$pk.loc-ipkx$pk.loc)),]
      allre[[i]]=cbind(re,Sid=isid,tick.max=max(v0[re[1,2]:re[1,3]],na.rm=TRUE),PkCl=currpk)
      next
    }
    
    re=.MGsimpleIntegr(x=lx,y=v,noise.local = rep(minNoise,length(lx)) ,snr.thresh = 1.01,span=ceiling(nspan/4)*2+1,minNoise = minNoise*1.01,v2alim = 0.8)
    if(nrow(re)){
      re=.MGinfctmergepk(re)
      re=re[which.min(abs(re$pk.loc-ipkx$pk.loc)),]
      allre[[i]]=cbind(re,Sid=isid,tick.max=max(v0[re[1,2]:re[1,3]],na.rm=TRUE),PkCl=currpk)
      next
    }
    
    #iloc=which.max(v0[lx])
    l=which(v>(max(v)/20) & v>quantile(v,.2))
    llpks=.GRsplist(l,l,d=1.1)
    llpkmax=suppressWarnings(sapply(llpks,function(x) max(v0[lx[x]],na.rm=T)))
    if(!all(is.infinite(llpkmax))){
      llpks=llpks[[which.max(llpkmax)]]
      llpkslx=lx[llpks]
      allre[[i]]=data.frame(pk.loc=llpkslx[which.max(v0[llpkslx])],
                            pk.left=min(llpkslx),pk.right=max(llpkslx),
                            pk.snr=max(v[llpks]),pk.int=max(v[llpks]),
                            pk.span=diff(range(llpks)),pk.cl=1,
                            Sid=isid,tick.max=max(v0[llpkslx],na.rm=T),PkCl=currpk)
      next
    }
    if(doPlot) cat("Missing in ",isid," / ",currpk," in ",main,"\n",sep="")
  }
  if(any(!sapply(allre,is.null))){
    addpks=do.call("rbind",allre)
    names(addpks)=gsub("^pk.","tick.",names(addpks))
    addpks=.MGinEICintrfct(addpks,m,m0,mbsl,mmz,parDeco,nspan)[[2]]
    addpks$PkCl=pks$PkCl[addpks$Pk]
    addpks$InDeco=0 ## for post processing analysis
    matpks=rbind(matpks,addpks)
  }
  
  
  ######################## Finalize ######################## 
  ### Set the rt back 
  sc2rt=as.numeric(rownames(m))*parDeco$psdrt
  rownames(matpks)=rownames(pks)=NULL
  for(i in names(matpks)[grep("^rt",names(matpks))]) matpks[,i]=round(sc2rt[matpks[,i]],3)
  for(i in c( "pk.loc", "pk.left" ,"pk.right")) pks[,i]=round(sc2rt[pks[,i]],3)
  if(nrow(apksrm)>0) for(i in c( "tick.loc", "tick.left" ,"tick.right")) apksrm[,i]=round(sc2rt[apksrm[,i]],3)
  matpks=matpks[order(matpks$Pk,matpks$Sid),]
  
  #### get the mz stats
  tmp=matpks[matpks$InDeco==1,]
  tmp=do.call("rbind",tapply(1:nrow(tmp),tmp$Pk,function(x){
    data.frame(mzmed=suppressWarnings(round(matrixStats:::weightedMedian(tmp$mzmed,tmp$int.sm,na.rm=T),6)),
               mzmin=round(min(tmp$mzmed),6),mzmax=round(max(tmp$mzmed),6))},simplify = FALSE))
  pks=data.frame(pks,tmp[match(pks$Pk,rownames(tmp)),])
  
  
  
  
  ############## Plot the whole stuff
  if(doPlot) .infctintegrplot(m,m0,parDeco,matpks,fac=10^6,rmz=range(mmz,na.rm=T),cols =colSa,sidsl = letSa,typs=typSa,main=main,v=renmfv,mpks=pks)
  
  ## return master peak??
  invisible(list(Pks=pks,MatPks=matpks,RmPks=apksrm))
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Internal function for formatting the final peak table
#' 
#' @param apks DF of peaks
#' @param m  transformed intensities matrix
#' @param m0 original intensities matrix
#' @param m0 baseline matrix
#' @param m0 original mz matrix
#' @param parDeco list of deconvolution paramters
#' @param nspan nspan
#' @keywords internal
#' 
#' @export
.MGinEICintrfct<-function(apks,m,m0,mbsl,mmz,parDeco,nspan){
  
  llpks=tapply(1:nrow(apks),paste0(apks$PkCl,";;",apks$Sid),list)
  
  sc2rt=as.numeric(rownames(m))*parDeco$psdrt
  minpk=matrix(FALSE,ncol=ncol(m),nrow=nrow(m),dimnames = dimnames(m))
  matpks=list()
  for(idx in names(llpks)){
    i=llpks[[idx]]
    lx=range(apks[i,1:3])
    lx=lx[1]:lx[2]
    lrt=sc2rt[lx]
    minpk[lx,apks$Sid[i]]=T
    v=m[lx,apks$Sid[i]]
    vb=mbsl[lx,apks$Sid[i]]
    # vvb=v-vb;vvb[vvb<0]=0
    v0=v00=m0[lx,apks$Sid[i]]
    v00[is.na(v0)]=parDeco$minNoiseMS1
    vmz=mmz[lx,apks$Sid[i]]
    if(all(is.na(vmz))) next
    lnna=which(!is.na(v0))
    lnna=lnna[order(abs(lnna-which.max(v)),-v0[lnna])]
    if(length(lnna)>nspan) lnna=lnna[1:nspan] 
    iap=lnna[which.max(v0[lnna])]
    
    matpks[[idx]]=data.frame(Sid=apks$Sid[i][1],Pk=apks$PkCl[i][1],
                             rtmin=min(lx),
                             rtmax=max(lx),
                             rt=lx[which.max(v)],
                             rtap=lx[iap],
                             npts=sum(!is.na(v0)),
                             mz=round(vmz[iap],6),
                             mzmed=suppressWarnings(round(matrixStats:::weightedMedian(vmz,v0,na.rm=T),6)),
                             snr.sm=round(max(v/vb),3),
                             int.sm=round(max(v),3),
                             int.ap=v0[iap],
                             area.sm=round(unname(.GRgetArea(lrt,v)[1]),3),
                             area.bsl=round(unname(.GRgetArea(lrt,vb)[1]),3),
                             area=round(.GRgetArea(lrt,v00-parDeco$minNoiseMS1)[1],3))
  }
  
  matpks=do.call("rbind",matpks)
  return(list(minpk,matpks))
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Internal function for collapsing peak table
#' 
#' @param tmp DF of peaks
#' @keywords internal
#' 
#' @export
.MGinfctmergepk=function(tmp)
  do.call("rbind",tapply(1:nrow(tmp),tmp$pk.cl,function(x){
    iap=which.max(tmp$pk.int[x])
    data.frame(pk.loc=tmp$pk.loc[iap], "pk.left"=min(tmp$pk.left[x]), "pk.right"=max(tmp$pk.right[x]),
               "pk.snr"=tmp$pk.snr[iap], "pk.int"=tmp$pk.int[iap],
               "pk.span"=diff(range(tmp$pk.left[x],tmp$pk.right[x])), "pk.cl"=tmp$pk.cl[x][1])
  }))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Internal function for baseline calculation
#' 
#' @param y  vector of intensities
#' @param span span
#' @param minNoise noise intensity level
#' @param n2pad left/right padding
#' @keywords internal
#' 
#' @export
.MGinfctonebsl<-function(y,span=11,minNoise=5000,n2pad=131){
  
  n0=length(y)
  y=c(rep(minNoise,n2pad-span*3-1),rep(y[1],span*3+1),y,rep(rev(y)[1],span*3+1),rep(minNoise,n2pad-span*3-1))
  bsl=GRMeta:::.GRbslrf(1:length(y),y,NoXP = NULL)
  bsl$fit[bsl$fit<minNoise]=minNoise
  bslscore <- (y - bsl$fit)/max(bsl$sigma, 10^-3)
  bslscore[which(abs(bslscore) > 10)] = sign(bslscore[which(abs(bslscore) > 10)]) * 10
  
  list(bsl$fit[n2pad+(1:n0)],bslscore[n2pad+(1:n0)])
}
