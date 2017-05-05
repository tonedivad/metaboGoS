### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name .MGmultiEICintgrPar
#' @title multiple integration of several eic
#' 
#' Multicore wrapper for .MGmultiEICintgr
#'
#' @param alleicmat eic: list of mz, y, and y2
#' @param parDeco deconvolution parameter list
#' @param l2excl vector of sample ids to be excluded when refining peak allocations
#' @param doPlot plotting at each step, vector of 0,1,2 - disabled if run in parallele
#' @param nSlaves num slaves
#' @import matrixStats
#' @return peak matrix data.frame
#' @export
.MGmultiEICintgrPar<-function(alleicmat,parDeco,l2excl=NULL,doPlot=1,nSlaves=1){
  
llre=list()
ll=names(alleicmat)#[111:120]
lperc=""


ll=ll[which(sapply(alleicmat[ll],function(x) max(x$y,na.rm=T))>=parDeco$minHeightMS1)]
if(length(ll)>20) lperc=ll[round(seq(1,length(ll),length=12)[2:11])]

if(nSlaves>1)   nSlaves=max(1, min(nSlaves,detectCores()-1))
if(nSlaves>1){
  clProc<-makeCluster(nSlaves)
  doParallel::registerDoParallel(clProc)
  cat(" -- registering ",nSlaves," clusters\n",sep="")
}

### Parallele bit
if(nSlaves>1)
   llre=foreach(idx = ll,.packages = c("metaboGoS"), .verbose =FALSE)  %dopar%{
    re=.MGmultiEICintgr(eic=alleicmat[[idx]],parDeco,doPlot = 0,l2excl = l2excl)
    if(is.null(re)) return(list())
    data.frame(RoiId=idx,re)
  }
## Serial bit
if(nSlaves<=1) for(idx in ll){
  re=.MGmultiEICintgr(eic=alleicmat[[idx]],parDeco,doPlot = doPlot,l2excl = l2excl)
  if(is.null(re)) next
  llre[[idx]]=data.frame(RoiId=idx,re)
}

if(nSlaves>1) stopCluster(clProc)

## combine results
llre=llre[sapply(llre,length)>0]
if(length(llre)==0){
  cat("Something went wrong no acceptable ROI were found!!!")
  return(obj)
}
ares=do.call("rbind",llre)
return(ares)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name .MGmultiEICintgr
#' @title multiple integration of several eic
#' 
#' stuff
#'
#' @param eic eic: list of mz, y, and y2
#' @param parDeco deconvolution parameter list
#' @param l2excl vector of sample ids to be excluded when refining peak allocations
#' @param doPlot plotting at each step, vector of 0,1,2
#' @import matrixStats
#' @return peak matrix data.frame
#' @export
.MGmultiEICintgr<-function(eic,parDeco=list(psdrt=0.00446,span=7,
                                            minHeightMS1=100000,minNoiseMS1=5000),l2excl=NULL,doPlot=0){

  m=t(eic$y2)#/1000000
  m0=t(eic$y)#/1000000
  mmz=t(eic$mz)#/1000000
  if(all(colnames(m)%in%l2excl)) return(NULL) ## if all to be excluded then stop
  maxm0=max(m0[,!colnames(m0)%in%l2excl],na.rm=T)
  if(maxm0<parDeco$minHeightMS1) return(NULL)
  span=parDeco$span
  minNoise=parDeco$minNoiseMS1
  
  nspan=floor(span/2)*2+1
  
  
bsllist=list(bsl=m,bslsc=m)
for(isid in colnames(m)){
  y=m[,isid]
  n2pad=131
  y=c(rep(minNoise,n2pad-span*3-1),rep(y[1],span*3+1),y,rep(rev(y)[1],span*3+1),rep(minNoise,n2pad-span*3-1))
  bsl=GRMeta:::.GRbslrf(1:length(y),y,NoXP = NULL)
  bsl$fit[bsl$fit<minNoise]=minNoise
  bslscore <- (y - bsl$fit)/max(bsl$sigma, 10^-3)
  bslscore[which(abs(bslscore) > 10)] = sign(bslscore[which(abs(bslscore) > 10)]) * 10
  
  bsllist$bsl[,isid] = bsl$fit[n2pad+(1:nrow(m))]
  bsllist$bslsc[,isid] = bslscore[n2pad+(1:nrow(m))]
  
}


######################## get the master peak list ######################## 
##     pks=.MGsimpleIntegr2(x=newx[,"rt"],y=newx[,"y2"],bsl = newx$bsl,bslscore = newx$bslc,snr.thresh = parDeco$sbslr,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8)

apks=list()
for(isid in colnames(m)){
  # x=1:nrow(m);y=m[,isid];bsl = bsllist$bsl[,isid];bslscore = bsllist$bslsc[,isid];snr.thresh =parDeco$sbsl;  span=nspan;minNoise = minNoise*1.01;v2alim = 0.8
  
  re=.MGsimpleIntegr2(x=1:nrow(m),y=m[,isid],bsl = bsllist$bsl[,isid],bslscore = bsllist$bslsc[,isid],snr.thresh =parDeco$sbslr,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8)
  if(nrow(re)==0) next
  re$Sid=isid
  apks[[isid]]=re
}
apks=do.call("rbind",apks)
if(is.null(apks)) return(NULL)
names(apks)=gsub("pk\\.","tick.",names(apks))
apks$PkCl=NA

##### Split large peak lists into chunks
llpk=.GRsplist(apks$tick.loc,d=6.1*nspan)
if(any(sapply(llpk,is.list))) llpk=unlist(llpk,recursive = F)
for(i in 1:length(llpk)) apks$PkCl[llpk[[i]]]=i

## only keep the ones with large 
l2k=as.numeric(names(which(tapply(apks$tick.int*2,apks$Pk,max)>parDeco$minHeightMS1 & tapply(apks$tick.snr,apks$Pk,max)>2)))
l2k=l2k[l2k%in%unique(apks$PkCl[!apks$Sid%in%l2excl])]
apks=apks[apks$Pk%in%l2k,]
if(nrow(apks)==0) return(NULL)
apks$PkCl=as.numeric(factor(apks$PkCl,names(sort(tapply(apks$tick.loc,apks$PkCl,mean)))))
apks=apks[order(apks$Sid,apks$PkCl),]
apks$PkCl2=1

## set of peaks that have redundacies/convolved
l2chk=names(which(tapply(apks$tick.loc,apks$PkCl,function(x) diff(range(x)))>nspan |
                    tapply(apks$Sid,apks$PkCl,function(x) any(duplicated(x)))))
i2chk=l2chk[1]

for(i2chk in l2chk){
  l2verif=which(apks$PkCl==i2chk)
  apks$PkCl2[l2verif]=.infctverif(apks[l2verif,],m=m,nspan = nspan,l2excl = l2excl,doplot=(3%in%doPlot))
}
newcl=paste(apks$PkCl,apks$PkCl2,sep=";")
apks$PkCl2=as.numeric(factor(newcl,names(sort(tapply(apks$tick.loc,newcl,mean,na.rm=T)))))

## diaply final alignment?
if(2%in%doPlot){
  cols=rep(brewer.pal(8,"Dark2"),length(unique(apks$Sid))/8)
  cols=cols[as.numeric(factor(apks$Sid,unique(apks$Sid)))]
  plot(apks$tick.loc,col=cols,pch=21,cex=.5)
abline(h=tapply(apks$tick.loc,apks$PkCl2,mean),lty=2,col="grey")
text(1:nrow(apks),apks$tick.loc,apks$PkCl2,col=cols)
}

######################## Form the final peak table ######################## 
llpks=tapply(1:nrow(apks),paste0(apks$PkCl2,";;",apks$Sid),list)

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
  vb=bsllist$bsl[lx,apks$Sid[i]]
  vvb=v-vb;vvb[vvb<0]=0
  v0=m0[lx,apks$Sid[i]]
  vmz=mmz[lx,apks$Sid[i]]
  lnna=which(!is.na(v0))
  lnna=lnna[order(abs(lnna-which.max(v)),-v0[lnna])]
  if(length(lnna)>nspan) lnna=lnna[1:nspan] 
  iap=lnna[which.max(v0[lnna])]
  
  matpks[[idx]]=data.frame(Sid=apks$Sid[i][1],Pk=apks$PkCl2[i][1],
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
                           area.sm=round(unname(.GRgetArea(lrt,vvb)[1]),3),
                           area=round(.GRgetArea(lrt[!is.na(v0)],v0[!is.na(v0)])[1],3))
}

matpks=do.call("rbind",matpks)
lupk=which(tapply(matpks$int.ap,matpks$Pk,max)>=parDeco$minHeightMS1 &
             tapply(matpks$snr.sm,matpks$Pk,max)>=2 &
             tapply(matpks$npts,matpks$Pk,max)>=floor(nspan/2))
if(length(lupk)==0) return(NULL)
matpks=matpks[matpks$Pk%in%names(lupk),]
lupk=sort(tapply(matpks$rt,matpks$Pk,median))
matpks$Pk=as.numeric(factor(matpks$Pk,names(lupk)))
matpks=matpks[order(matpks$Pk,matpks$Sid),]
matpks$InDeco=1 ## for post processing analysis

######################## Fix missing ######################## 
## these are due to peaks below snr
pksum=cbind(round(tapply(matpks$rt,matpks$Pk,quantile,0.5)),
            floor(tapply(matpks$rtmin,matpks$Pk,quantile,0.25)),
            ceiling(tapply(matpks$rtmax,matpks$Pk,quantile,0.75)))
lmiss=names(which(table(matpks$Sid)!=max(matpks$Pk)))
pks2add=list()
for(isid in lmiss){
  lin=unique(matpks$Pk)[!unique(matpks$Pk)%in%matpks$Pk[matpks$Sid==isid]]
  for(ipk in lin){
    lx=pksum[ipk,2]:pksum[ipk,3]
    v=m[lx,isid]
    iv=minpk[lx,isid]
    l=which(matpks$Pk>ipk & matpks$Sid==isid);if(length(l)) iv[lx>min(matpks$rt[l])]=TRUE
    l=which(matpks$Pk<ipk & matpks$Sid==isid);if(length(l)) iv[lx<max(matpks$rt[l])]=TRUE
    if(all(v[!iv]<(minNoise*2))) next
    
    # x=lx;y=v;noise.local = rep(minNoise,length(lx)) ;snr.thresh = 1.01;span=ceiling(nspan/2);minNoise = minNoise*1.01;v2alim = 0.8
    
    re=.MGsimpleIntegr(x=lx,y=v,noise.local = rep(minNoise,length(lx)) ,snr.thresh = 1.01,span=ceiling(nspan/2),minNoise = minNoise*1.01,v2alim = 0.8)
    if(nrow(re)==0) next
    ##  make sure the new peak is in no peak area
    re=re[which(re$pk.loc%in%lx[!iv]),]
    if(nrow(re)==0) next
    re=re[which.min(abs(re$pk.loc-pksum[ipk,1])),]
    v0=m0[re[,2]:re[,3],isid]
    if(all(is.na(v0))) next
    vmz0=mmz[re[,2]:re[,3],isid]
    ap0=which.max(v0)
    vvb=(m[,isid]-bsllist$bsl)[re[,2]:re[,3]];vvb[vvb<0]=0
    pks2add=c(pks2add,list(data.frame(Sid=isid,Pk=ipk,
               rtmin=re[,2],rtmax=re[,3],rt=re[,1],rtap=ap0+re[,2]-1,
               npts=sum(!is.na(v0)),
    mz=round(vmz0[ap0],6),mzmed=round(matrixStats:::weightedMedian(vmz0,v0,na.rm=T),6),
    snr.sm=re$pk.snr,int.sm=re$pk.int,int.ap=v0[ap0],
    area.sm=round(unname(.GRgetArea(sc2rt[re[,2]:re[,3]],vvb)[1]),3),
    area=round(.GRgetArea(sc2rt[re[,2]:re[,3]][!is.na(v0)],v0[!is.na(v0)])[1],3),InDeco=0)))
  }
}
pks2add=do.call("rbind",pks2add)
if(!is.null(pks2add)) matpks=rbind(matpks,pks2add)


######################## Finalize ######################## 
### Set the rt back 
rownames(matpks)=NULL
for(i in names(matpks)[grep("^rt",names(matpks))]) matpks[,i]=round(sc2rt[matpks[,i]],3)
matpks=matpks[order(matpks$Pk,matpks$Sid),]

if(1%in%doPlot) .infctintegrplot(m,m0,parDeco,matpks,fac=10^6,range(mmz,na.rm=T))

############## Plot the whole stuff

invisible(matpks)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Internal function for plotting the outcome multiple integration of EICs
#' 
#' @param ipks DF of peaks
#' @param m EIC matrix
#' @param m0 EIC matrix using the original datapoints
#' @param parDeco list of deconvolution paramters
#' @param matpks matrix of integrated peaks
#' @param fac intensity factor to make intensity easier to read
#' @keywords internal
#' 
#' @export
.infctintegrplot<-function(m,m0,parDeco,matpks,fac=10^6,rmz){
  
  lso=order(-colSums(m))
  m=sqrt(m[,lso,drop=F]/fac)
  m0=sqrt(m0[,lso,drop=F]/fac)
  lrt=as.numeric(rownames(m))*parDeco$psdrt
  sids=colnames(m)
  sidsl=rep(c(letters,LETTERS),ceiling(length(sids)/52))[1:length(sids)]
  names(sidsl)=sids
  
  xl=pretty(lrt)
  ylim=pretty(c(0,max(m,na.rm=T),max(m0,na.rm=T),sqrt(matpks$int.sm/fac)))
  cols=rep(brewer.pal(8,"Dark2"),ceiling(ncol(m)/8))[1:ncol(m)]
  names(cols)=colnames(m)
  
  #par(mfrow=c(2,1),mar=c(5,4,1,.1))
  par(mfrow=c(2,1),mar=c(1,4,1,.1))
  matplot(lrt,m,typ="l",ylim=range(ylim),xlim=range(xl),#main=iroi,
          col=cols,axes=F,xlab="rt",ylab="Sqrt(Int)")
  # matplot(lrt,m0,typ="p",ylim=ylim,col=metaData[sids,]$Cols,pch=16,add=T)
  axis(2,at=ylim,ylim^2,las=2,pos=min(xl))
 # axis(1,at=xl,pos=min(ylim))
  abline(v=xl,lty=2,col="grey")
  legs=paste0(sids," (",colSums(!is.na(m0)),")")
  legend("topright",legs,col=cols,cex=.5,bty="n",ncol=1,pch=sidsl,pt.cex=0.5)
  cols=rep(brewer.pal(8,"Dark2"),ceiling(max(matpks$Pk)/8))
  text(matpks$rt,sqrt(matpks$int.sm/fac),sidsl[matpks$Sid],col=cols[matpks$Pk])
  abline(v=tapply(matpks$rtmax,matpks$Pk,quantile,.75)-parDeco$psdrt*2,col=brewer.pal(8,"Dark2")[1:max(matpks$Pk)],lty=2,lwd=2)
  abline(v=tapply(matpks$rtmin,matpks$Pk,quantile,.25)+parDeco$psdrt*2,col=brewer.pal(8,"Dark2")[1:max(matpks$Pk)],lty=3,lwd=2)
  
  ylim=pretty(rmz)
  par(mar=c(5,4,0,.1))
  plot(range(matpks$rt),range(matpks$mz),ylim=range(ylim),xlim=range(xl),axes=F,xlab="rt",ylab="",col=cols[matpks$Pk],pch=16,cex=0) #,main=iroi
  text(matpks$rt,matpks$mz,sidsl[matpks$Sid],col=cols[matpks$Pk])
  
  abline(v=tapply(matpks$rtmax,matpks$Pk,quantile,.75)-parDeco$psdrt*2,col=brewer.pal(8,"Dark2")[1:max(matpks$Pk)],lty=2,lwd=2)
  abline(v=tapply(matpks$rtmin,matpks$Pk,quantile,.25)+parDeco$psdrt*2,col=brewer.pal(8,"Dark2")[1:max(matpks$Pk)],lty=3,lwd=2)
  abline(h=ylim,lty=2,col="grey")
  abline(v=xl,lty=2,col="grey")
  axis(2,at=ylim,las=2,pos=xl[1]);axis(1,at=xl,pos=ylim[1])
  par(mfrow=c(1,1),mar=c(5,4,1,.1))
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' Internal function for multiple integration of EICs
#' 
#' @param ipks DF of peaks
#' @param m EIC matrix
#' @param nspan span as in scan
#' @param l2excl vector of sample ids to be excluded from the refinmet
#' @param doplot TRUE/FALSE plot the outcome of the refinment process
#' @import NMF
#' @keywords internal
#' 
#' @export
.infctverif<-function(ipks,m,nspan,l2excl=NULL,doplot=FALSE){
  lsc=range(ipks[,1:3])+c(-1,1)*nspan
  lsc[1]=max(1,lsc[1]);lsc[2]=min(nrow(m),lsc[2]);lsc=lsc[1]:lsc[2]
  
  lusamp=unique(ipks$Sid)
  lusamp=lusamp[!lusamp%in%l2excl]
  if(length(lusamp)>1) re=NMF:::loadings(NMF:::nmf(m[lsc,lusamp,drop=F], 1, 'snmf/l',seed='ica')) else re=m[lsc,lusamp,drop=F]
  #tabres=.GRcompSeg(m[lsc,isamp],bsllist$bsl[lsc,isamp],list(nspan=nspan),5000)
  
  linspan=rev(seq(3,nspan*2+1,2))
  inspan=1
  mpks=.GRmsPeakSimple(lsc,re[,1],span = linspan[inspan])
  while(nrow(mpks)==0 & inspan<length(linspan)){
    inspan=inspan+1
    mpks=.GRmsPeakSimple(lsc,re[,1],span = linspan[inspan])
  }
  if(nrow(mpks)==0) return(ipks$PkCl2) ## pb with some peak going down???
  linmpks=lapply(1:nrow(mpks),function(x) mpks[x,'mass.left']:mpks[x,'mass.right'])
  
  mperc=do.call("rbind",lapply(1:nrow(ipks),function(i){
    lsci=ipks[i,2]:ipks[i,3]
    sapply(linmpks,function(y) sum(m[lsci[lsci%in%y],ipks$Sid[i]]))/sum(m[lsci,ipks$Sid[i]])
  }))
  v=apply(mperc,1,which.max)
  v[apply(mperc,1,max)<0.1]=NA
  if(!doplot) return(v)
  
  plot(ipks$tick.loc,col=brewer.pal(11,"Spectral")[as.numeric(cut(apply(mperc,1,max),seq(0,1,.1)))],
       pch=15+apply(mperc,1,which.max),ylim=range(lsc))
  text(1:nrow(ipks),ipks$tick.loc,apply(mperc,1,which.max))
  abline(h=tapply(ipks$tick.loc,ipks$PkCl,mean))
  abline(h=mpks$mass.loc,col=2)
  
  invisible(v)
}

#############################################################################
##
#' Gather all eics from multiple files
#' @param lfiles list of files containing the EIC : sholud allxeic
#' @param eicmat data.frame containing EI?C roi infos
#' @param outfile file where alleicmat/eicmat/lfiles could be stored
#' @return list of aligned eics
#' @export
.MGgatherEICfromFiles<-function(lfiles,eicmat,outfile=NULL){
  
  lfiles=lfiles[file.exists(lfiles)]
  if(length(lfiles)==0) stop('No files found')
  if(is.null(names(lfiles))) names(lfiles)=paste0("S",1:length(lfiles))
  cat(" ++ found ",length(lfiles)," sample files\n",sep="")
  alleic=list()
  for(isid in names(lfiles)){
    cat(isid,"- ")
    load(lfiles[isid])
    allxeic=allxeic[names(allxeic)%in%eicmat$RoiId]
    alleic[[isid]]=allxeic
  }
  cat("\n")
  llre=list()
  ll=eicmat$RoiId#[111:120]
  lperc=ll[round(seq(1,length(ll),length=12)[2:11])]
  
  cat(" ++",nrow(eicmat)," EICs to merge:\n",sep="")
  
  alleicmat=list()
  eicmat$scmax=eicmat$scmin=eicmat$nrtmax=eicmat$nrtmin=NA
  for(iroi in ll){
    if(iroi %in% lperc) cat(iroi,"(",which(ll==iroi),") ",sep="")
    ieic=lapply(alleic,function(x) x[[iroi]])
    ieic=ieic[!sapply(ieic,is.null)]
    if(length(ieic)==0) next
    lsc=max(sapply(ieic,function(x) min(x[,"nscan"]))):min(sapply(ieic,function(x) max(x[,"nscan"])))
    eicmat[which(eicmat$RoiId==iroi),c("scmin","scmax")]=range(lsc)
    eicmat[which(eicmat$RoiId==iroi),c("nrtmin","nrtmax")]=c(max(sapply(ieic,function(x) min(x[,"rt"]))),min(sapply(ieic,function(x) max(x[,"rt"]))))
    ieic=lapply(ieic,function(x) x[match(lsc,x[,"nscan"]),])
    ieic=list(y=t(sapply(ieic,function(x) x[,"y"])),
              mz=t(sapply(ieic,function(x) x[,"mz"])),
              y2=t(sapply(ieic,function(x) x[,"y2"])))
    alleicmat[[iroi]]=lapply(ieic,function(x){colnames(x)=lsc;x})
    
  }
  if(!is.null(outfile)){
    cat(" ++ saving to ",outfile,"\n",sep="")
    save(file=outfile,alleicmat,eicmat,lfiles)
  }
  
  invisible(alleicmat)
}

