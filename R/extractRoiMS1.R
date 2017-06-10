### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name extractRoiMS1
#' @title extractRoiMS1
#' 
#' Refine/integra ROI 
#'
#' @param samp xml file or xcms object
#' @param parDeco deconvolution parameter list
#' @param blList list of blank filles xcmsRaw object 
#' @param lmzExcl list of m/z excluded form the blank correction
#' @param doPlot plotting at each step
#' @param mzrange mzrange
#' @param useMinNoise useMinNoise
#' @return RoiInfos / PeakInfos data.frame
#' @export
extractRoiMS1<-function(samp,parDeco,blList=NULL,lmzExcl=c(),stepStop=3,doPlot=TRUE,mzrange=c(-Inf,Inf),useMinNoise=TRUE){
  
  # blList=NULL;lmzExcl=NULL;stepStop=3;mzrange=c(-Inf,Inf);doPlot=TRUE
  

## samp is an xcms file
if("xcmsRaw"%in%class(samp)){
  xRaw=samp
  samp=unclass(xRaw@filepath)
  cat(" As provided ",samp,"\n",sep="")
}else{
  cat(" Parsing ",samp,"\n",sep="")
  xRaw <- xcmsRaw(samp, profmethod = "bin", profparam = list(step = 0.0001), profstep = 0)
}

sctime=(xRaw@scantime)/60
rangeSc=c(1,1,length(sctime),length(sctime))
if(!all(is.na(parDeco$rtlim))){
  if(length(parDeco$rtlim)==2) rangeSc[2:3]=range(which(sctime>=parDeco$rtlim[1] & sctime<=parDeco$rtlim[2]))
  if(length(parDeco$rtlim)==4){
    rangeSc[2:3]=range(which(sctime>=parDeco$rtlim[2] & sctime<=parDeco$rtlim[3]))
    rangeSc[c(1,4)]=range(which(sctime>=parDeco$rtlim[1] & sctime<=parDeco$rtlim[4]))
  }
}



######## Step 1 ######## ######## 
## looking for ROIs with consecutive
strt=Sys.time()
cat(" Looking for ROIs",sep="")
ncons=max(3,parDeco$ncons-1)
nhei=parDeco$minNoiseMS1*parDeco$sbr
if(!useMinNoise & !is.null(parDeco$minHeightMS1cons)) nhei=parDeco$minHeightMS1cons

ROI.list=GRMeta:::.GRfindROIs(xRaw, parDeco$ppm,prefilter = c(ncons,nhei),maxskip =0)
# ROI.list <- xcms:::findmzROI(xRaw, dev = parDeco$ppm * 1e-06, minCentroids = 0,prefilter = c(max(2,parDeco$ncons-1),parDeco$minNoise*parDeco$sbr), noise = parDeco$minZero)

ROImat=do.call("rbind",lapply(ROI.list,unlist))
ROImat=ROImat[order(ROImat[,1]),]
l=which(ROImat[,"scmax"]>=rangeSc[1] & ROImat[,"scmin"]<=rangeSc[4])
ROImat=ROImat[l,,drop=F]
ROImat[,"intensity"]=apply(ROImat,1,function(x) max(rawEIC(xRaw,mzrange=x[2:3],scanrange=x[4:5])$intensity))
cat(" 0 --> found ",nrow(ROImat)," ROIs [+",round(Sys.time()-strt,2),"sec.]\n",sep="")
if(stepStop==0) return(ROImat)


######## Step 2 ######## ######## 
## merge close m/z and potentially overlapping peaks
strt=Sys.time()
llppm<-.GRsplistMZ(ROImat[,1],dppm = parDeco$ppm,dmz=parDeco$dmz,typ="min") ## min distance b/w dppm/dmz
#table(sapply(llppm,length))
ROImat1<-do.call("rbind",lapply(llppm,function(x)
  c("mz"=median(ROImat[x,"mz"]),"mzmin"=min(ROImat[x,"mzmin"]),"mzmax"=max(ROImat[x,"mzmax"]),
    "scmin"=min(ROImat[x,"scmin"]),"scmax"=max(ROImat[x,"scmax"]),"length"=0,intensity=max(ROImat[x,"intensity"]))))
l=which(ROImat1[,"scmax"]>=rangeSc[2] & ROImat1[,"scmin"]<=rangeSc[3] & ROImat1[,"mzmax"]>=(mzrange[1]-0.3) & ROImat1[,"mzmin"]<=(mzrange[2]+0.3))
ROImat1=ROImat1[l,,drop=F]
cat(" 1 --> reduced to ",nrow(ROImat1)," ROIs when merging close m/z in {",rangeSc[2],"-",rangeSc[3],"} (",
    round(mzrange[1]-.3,4),";",round(mzrange[2]+.3,4),") [+",round(Sys.time()-strt,2),"sec.]\n",sep="")
if(stepStop==1) return(ROImat1)

######## Step 3 ######## ######## 
### widen window to get the maximum sum intensity 
strt=Sys.time()
ROImat2=lapply(1:nrow(ROImat1),function(ix){
  xroi=ROImat1[ix,]
  lmz=range(c(xroi[c(1,1,2,3)]*(1+1*c(-2,2,-1,1)*parDeco$ppm*10^-6),xroi[2:3]+c(-1,1)*parDeco$dmz))+c(-1.1,1.1)*parDeco$dmz
  lrt=sctime[xroi[4:5]]+c(-1,1)*parDeco$drt+ncons*parDeco$psdrt
  lrt[1]=max(parDeco$rtlim[1],lrt[1])*60
  lrt[2]=min(parDeco$rtlim[4],lrt[2])*60
  xr=rawEIC(xRaw, mzrange =lmz,rtrange = lrt)
  # xr$scan=sctime[xr$scan]
  c(xroi[1],lmz[1],lmz[2],range(xr$scan),0,max(xr$intensity),xr$scan[which.max(xr$intensity)])
})
ROImat2=do.call("rbind",ROImat2)
colnames(ROImat2)=c(colnames(ROImat1),"scap")

l2k=which(ROImat2[,"intensity"]>=parDeco$minHeightMS1cons)
ROImat2=ROImat2[l2k,]
l=which(ROImat2[,"scmax"]>=rangeSc[2] & ROImat2[,"scmin"]<=rangeSc[3] & ROImat2[,"mzmax"]>=(mzrange[1]-0.1) & ROImat2[,"mzmin"]<=(mzrange[2]+0.1))
ROImat2=ROImat2[l,,drop=F]
ROImat2[,4:5]=.GRchkrange(ROImat2[,4:5,drop=FALSE],rangeSc = rangeSc[c(1,4)])
lso=order(round(rowMeans(ROImat2[,2:3,drop=FALSE]),4),rowMeans(ROImat2[,4:5,drop=FALSE]))
ROImat2=ROImat2[lso,,drop=FALSE]
rownames(ROImat2)=1:nrow(ROImat2)
cat(" 2 --> reduced to ",nrow(ROImat2)," ROIs above ",parDeco$minHeightMS1cons," cps in {",rangeSc[2],"-",rangeSc[3],
    "} (",round(mzrange[1]-.1,4),";",round(mzrange[2]+.1,4),") [+",round(Sys.time()-strt,2),"sec.]\n",sep="")
ROImat2=cbind(ROImat2,coda=1,roi=1:nrow(ROImat2))
if(stepStop==2) return(ROImat2)

######## Step 4 ######## ######## 
### overlap  and integrate
llover=.GRisover(ROImat2[,2],ROImat2[,3],retOne = T,thr = parDeco$dmz*2)
cat(" 3: refining/ingretating ",length(llover)," groups of rois ",sep="")

newrtx=(1:ceiling(max(xRaw@scantime/60/parDeco$psdrt)))*parDeco$psdrt
sc2nrt=apply(abs(outer(newrtx,xRaw@scantime/60,"-")),2,which.min)
sc2nrt=cbind(scan=1:length(xRaw@scantime),nscan=sc2nrt,nrt=round(newrtx[sc2nrt],5))

ROIlist=list()
addrt=parDeco$psdrt*(2*parDeco$span2+3)+parDeco$drt

for(ill in (1:length(llover))){
  xroi=ROImat2[llover[[ill]],,drop=F]
  if(ill%%100==0) cat(ill," ")
  
  lmz=range(xroi[,1:3])
  lsc=range(xroi[,4:5])
  eic=data.frame(.GRrawMat(xRaw,mzrange = lmz, scanrange = lsc,padsc =T,naVal = parDeco$minNoiseMS1),Sid="S")
  resref=.MGreframeOneROI(eic,parDeco,minHeight1 = parDeco$minHeightMS1cons,minHeight2 = parDeco$minHeightMS1,addrt,doPlot = doPlot,main = ill)
  if(is.null(resref)) next
  repks=list()
  for(ieic in 1:length(resref)){
    xeic=resref[[ieic]]
    xeic=cbind(xeic,sc2nrt[xeic[,"scan"],2:3])
    
    lmz=range(xeic$mz,na.rm = T)+c(-1,1)*parDeco$dmz
    lrt=range(xeic$rt,na.rm = T)+c(-2,2)*parDeco$psdrt
    eicbl=NULL
    if(length(blList)>0){
      eicbl=do.call("rbind",lapply(1:length(blList),function(ibl){
        i=cbind(GRMeta:::.GRrawMat(blList[[ibl]],mzrange = lmz, rtrange = lrt*60+c(-10,10),padsc =T),blid=ibl)
        if(sum(!is.na(i[,"mz"]))<5) return(NULL)
        i
      }))
    }
    isbrlim=ifelse(any((lmzExcl-lmz[1]-parDeco$dmz)>0 & (lmz[2]+parDeco$dmz-lmzExcl)>0),0,parDeco$sbr) ## if in excl
    repks[[ieic]]=.MGsingleIntegr(xeic,eicbl=eicbl,parDeco,isbrlim,doPlot=doPlot)
  }
  ROIlist[[ill]]=do.call("rbind",repks)
  
}
ROImat3=do.call("rbind",ROIlist)

cat("\n   --> reduced to ",nrow(ROImat3),
    "/",length(unique(ROImat3$cl.id)),
    "/",length(unique(ROImat3$roi.id))," peaks/peak clusters/ROIs after refining [+",round(Sys.time()-strt,2),"sec.]\n",sep="")

return(ROImat3)
}
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 



######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 

.MGsingleIntegr<-function(xeic,eicbl=NULL,parDeco,sbrlim,doPlot=T){
  
  nspan=parDeco$span
  span2=parDeco$span2
  
  newx=.MGfill(xeic,drt = parDeco$psdrt,span = parDeco$span,bw = parDeco$bw,minNoise = parDeco$minNoiseMS1,minpts = parDeco$ncons)
  
  
  ## add roi informations
  roiinf=data.frame(roi.id=sprintf("R%.4f@%.2f-%.2f",median(newx[,"mz"],na.rm=T),min(newx[,"rt"],na.rm=T),max(newx[,"rt"],na.rm=T)),
                    roi.coda=.GRcoda(newx$y2),
                    roi.mz50=unname(median(newx[,"mz"],na.rm=T)),
                    roi.mz10=unname(quantile(newx[,"mz"],.1,na.rm=T)),
                    roi.mz90=unname(quantile(newx[,"mz"],.9,na.rm=T)))
  imain=paste0(roiinf$roi.id," MCQ=",round(roiinf$roi.coda,3))
  
  
  rtr0=range(newx$rt)
  
  
  ##### check blank from and align to newx$x
  newx$ybl=rep(parDeco$minNoiseMS1,nrow(newx))
  if(!is.null(eicbl)){
    tmpbl=eicbl
    # l2excl=which(tmpbl[,"mz"]<mzr[1] | tmpbl[,"mz"]>mzr[2])
    # if(length(l2excl)) tmpbl[l2excl,c("mz","y")]=NA
    tmpbl=tmpbl[order(tmpbl[,"scan"],!is.na(tmpbl[,"mz"])),]
    tmpbl=tmpbl[which(!(duplicated(tmpbl[,"scan"]) & is.na(tmpbl[,"mz"]))),]
    tmpbl=tmpbl[tmpbl[,"blid"]%in%names(which(tapply(!is.na(tmpbl[,"y"]),tmpbl[,"blid"],sum)>10)),,drop=F]
    if(nrow(tmpbl)>1){
      tmpbl[which(is.na(tmpbl[,"y"]) | tmpbl[,"y"]<parDeco$minNoiseMS1),"y"]=parDeco$minNoiseMS1
      tmpbl<-try(do.call("cbind",tapply(1:nrow(tmpbl),tmpbl[,"blid"],function(xbl){
        ybl = .GRasysm(tmpbl[xbl,"y"], p = 0.5, lambda = 10^5)
        ybl[ybl <= parDeco$minNoiseMS1 | is.na(ybl)] = parDeco$minNoiseMS1
        ybl=approx(tmpbl[xbl,"rt"],ybl,newx$rt,yleft = ybl[1],yright = rev(ybl)[1])$y})),TRUE)
      if(!"try-error"%in%class(tmpbl)) newx$ybl=apply(tmpbl,1,median)
    }
  }
  
  pks=.MGsimpleIntegr2(x=newx$rt,y=newx$y2,bsl = newx$bsl,bslscore = newx$bsl,snr.thresh =parDeco$sbslr,span=nspan,
                       minNoise = parDeco$minNoiseMS1*1.01,v2alim = 0.8,span2 = span2+2)
  if(nrow(pks)==0){
    if(doPlot) .infctintegrplotOne(newx,parDeco,matpks = NULL,imain = imain)
    return(NULL)
  }
  
  ## compute the SBR around the apex
  pks$pk.sbr=round(sapply(pks$pk.loc,function(i) max((newx$y2/newx$ybl)[which(abs(newx$rt-i)<=(2*parDeco$psdrt))])),3)
  pks$pk.snr=round(sapply(pks$pk.loc,function(i) max((newx$y2/newx$bsl)[which.min(abs(newx$rt-i))])),3)
  
  #sbrlim=ifelse(any(lExcl>=irmz[1] & lExcl<=irmz[2]),0,parDeco$sbr)
  l2k=as.numeric(names(which(tapply(pks$pk.sbr,pks$pk.cl,max)>=sbrlim &
                               tapply(pks$pk.snr,pks$pk.cl,max)>=parDeco$sbslr &
                               tapply(pks$pk.int,pks$pk.cl,max)>parDeco$minHeightMS1)))
  pks=pks[pks$pk.cl%in%l2k,,drop=F]
  if(nrow(pks)==0){
    if(doPlot) .infctintegrplotOne(newx,parDeco,matpks = NULL,imain = imain)
    return(NULL)
  }
  
  pks$pk.cl=as.numeric(factor(pks$pk.cl))
  
  ## compute coda/roi stats
  pks$pk.roi=rep(1,nrow(pks))
  # if(nrow(pks)>1)
  #   for(i in 2:(nrow(pks))) pks$pk.roi[i]=pks$pk.roi[i-1]+ifelse((pks$pk.left[i]-pks$pk.right[i-1])>=parDeco$drt,1,0)
  
  ## integrate peaks
  newx$y2v=newx$y2-newx$bsl
  newx$y2v[newx$y2v<0]=0
  
  idx=1
  matpks=list()
  for(idx in 1:nrow(pks)){
    tmp=newx[which(newx$rt>=pks$pk.left[idx] & newx$rt<=pks$pk.right[idx]),]
    
    lnna=which(!is.na(tmp$y))
    if(length(lnna)==0) next
    lnna=lnna[order(abs(lnna-which.max(tmp$y2)),-tmp$y[lnna])]
    if(length(lnna)>nspan) lnna=lnna[1:nspan] 
    iap=lnna[which.max(tmp$y[lnna])]
    
    matpks[[idx]]=data.frame(id=idx,
                             pk.id=sprintf("P%.5f@%.3f",tmp$mz[iap],tmp$rt[iap]),
                             pk.rtmin=min(tmp$rt),
                             pk.rtmax=max(tmp$rt),
                             pk.rt=tmp$rt[which.max(tmp$y2)],
                             pk.rtap=tmp$rt[iap],
                             pk.int.ap=tmp$y[iap],
                             pk.int.apsm=max(tmp$y2),
                             pk.intmin=tmp$y2[1],
                             pk.bslmin=tmp$bsl[1],
                             pk.intmax=rev(tmp$y2)[1],
                             pk.bslmax=rev(tmp$bsl)[1],
                             pk.npts=sum(!is.na(tmp$y)),
                             pk.mz=round(tmp$mz[iap],6),
                             pk.mzmin=round(min(tmp$mz,na.rm = T),6),
                             pk.mzmax=round(max(tmp$mz,na.rm = T),6),
                             pk.mzmed=suppressWarnings(round(matrixStats:::weightedMedian(tmp$mz,tmp$y,na.rm=T),6)),
                             pk.area.sm=round(unname(.GRgetArea(tmp$rt,tmp$y2)[1]),3),
                             pk.area=round(.GRgetArea(tmp$rt[!is.na(tmp$y)],tmp$y[!is.na(tmp$y)])[1],3),
                             pk.snr=tmp$y[iap]/tmp$bsl[iap])
  }
  matpks=do.call("rbind",matpks)
  
  matpks=cbind(matpks,pk.sbr=pks[matpks$id,c("pk.sbr")],pks[matpks$id,grep("^roi.",names(pks))])
  matpks$id=paste0("R",pks$pk.roi,"C",pks$pk.cl)[matpks$id]
  matpks$cl.id=tapply(1:nrow(matpks),matpks$id,function(x) sprintf("C%.5f@%.3f",median(matpks$pk.mz[x]),median(matpks$pk.rtap[x])))[matpks$id]
  matpks=cbind(matpks,roiinf[rep(1,nrow(matpks)),])
  if(doPlot) .infctintegrplotOne(newx,parDeco,matpks = matpks,imain = imain)
  matpks
}
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 


######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
.infctintegrplotOne<-function(newx,parDeco,matpks=NULL,imain=NULL){
  
  yl=pretty(range(sqrt(c(0,newx$y2,newx$bsl*parDeco$sbslr,newx$y,newx$ybl*parDeco$sbr)),na.rm=T))
  # yl=pretty(sqrt(range(pretty(yl^2))))
  xl=pretty(newx$rt,n=11)
  par(mfrow=c(2,1))
  par(mar=c(1,5,2,0))
  if(is.null(imain)) imain=sprintf("[%.4f;%.4f] {%.2f-%.2f} MCQ=%.2f",min(newx[,"mz"],na.rm=T),max(newx[,"mz"],na.rm=T),
                                   min(newx[,"rt"],na.rm=T),max(newx[,"rt"],na.rm=T),.GRcoda(newx$y2))
  plot(newx[,"rt"],sqrt(newx[,"y"]),pch=16,main=imain,xlab="",ylab="Height*10^6",
       xlim=range(xl),ylim=range(yl),axes=F)
  axis(2,at=yl,yl^2/1000000,las=2,pos=min(xl))
  #  axis(1,at=xl,pos=min(yl))
  #      abline(v=c(unique(matpks$pk.rtmax),unique(matpks$pk.rtmin)),col="grey",lwd=2)
  lines(newx$rt,sqrt(newx$y2),col=2,lwd=2)
  lines(newx$rt,sqrt(newx$bsl),col=4,lwd=2,lty=1)
  lines(newx$rt,sqrt(newx$bsl*parDeco$sbslr),col=4,lwd=1,lty=4)
  lines(newx$rt,sqrt(newx$ybl*parDeco$sbr),col=3,lwd=1,lty=2)
  lines(newx$rt,sqrt(newx$ybl),col=3,lwd=2,lty=1)
  abline(h=sqrt(parDeco$minHeightMS1))
  abline(h=sqrt(parDeco$minHeightMS1cons),lty=3)
  
  if(!is.null(matpks)){
    cols2 = rep(brewer.pal(8, "Dark2"), ceiling(nrow(matpks)/7))
    points(matpks$pk.rtap,sqrt(matpks$pk.int.ap),pch=19,col=cols2,cex=2)
    abline(v=matpks$pk.rtmin+parDeco$psdrt/2,lty=3,col=cols2)
    abline(v=matpks$pk.rtmax-parDeco$psdrt/2,lty=4,col=cols2)
  }
  
  
  ################
  par(mar=c(4,5,0,0))
  yl=pretty(range(newx$mz,na.rm = T)*(1+c(-1,1)*2.1*10^-6))
  xl=pretty(newx$rt,n=11)
  
  cols=rep("black",nrow(newx))
  cols[which(newx[,"y"]<(parDeco$minHeightMS1))]="grey70"
  cols[which(newx[,"y"]<(parDeco$minHeightMS1cons))]="grey50"
  cols[which(newx[,"y"]<(newx$ybl*parDeco$sbr))]="grey30"
  
  plot(newx[,"rt"],newx[,"mz"],pch=16,main=NULL,xlab="Retention time",ylab="",xlim=range(xl),ylim=range(yl),axes=F,col=cols)
  axis(2,at=yl,las=2,pos=min(xl))
  nv=ceiling(diff(range(yl))/(median(yl)*10^-6))
  abline(h=(-nv:nv)*median(yl)*10^-6+median(yl),lty=3,col="grey")
  axis(1,at=xl,pos=min(yl))
  #      abline(v=c(unique(matpks$pk.rtmax),unique(matpks$pk.rtmin)),col="grey",lwd=2)
  if(!is.null(matpks)){
    points(matpks$pk.rtap,matpks$pk.mz,pch=19,col=cols2,cex=2)
    abline(v=matpks$pk.rtmin+parDeco$psdrt/2,lty=3,col=cols2)
    abline(v=matpks$pk.rtmax-parDeco$psdrt/2,lty=4,col=cols2)
  }
  par(mfrow=c(1,1),mar=c(5.1 ,4.1 ,4.1, 2.1))
}
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 


######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
.MGreframeOneROI<-function(eic,parDeco,minHeight1 = parDeco$minHeightMS1/2,minHeight2 = parDeco$minHeightMS1,
                           addrt=2*parDeco$psdrt*parDeco$span+parDeco$drt*2,doPlot=TRUE,main=""){
  #    for(ix in names(leics)){
  #     eic=leics[[ix]]
  
  if(!"Sid"%in%names(eic)) eic$Sid="s"
  l=which(!is.na(eic$mz))
  ieic=eic[l,]
  rmz=range(ieic[,"mz"])
  mindmz=min(c(mean(rmz)*parDeco$ppm*10^-6,parDeco$dmz))/2
  
  if(doPlot){
    ymz=range(pretty(range(eic$mz,na.rm=TRUE)+c(-0.001,0.001)))
    yrt=range(pretty(range(eic$rt,na.rm=TRUE)+c(-0.2,0.2)))
    l=which(!is.na(eic$mz))
    k <- MASS::kde2d(eic$mz[l], eic$rt[l],h = c(0.001,0.05),lims = c(ymz,yrt))
    #      k$z=k$z #*sum(!is.na(eic$mz))/sum(k$z)
    k$z[k$z<0.001]=NA
    image(k, col=colorRampPalette(rev(brewer.pal(9,'Greys')))(11)[11:6],main=main)
    abline(v=seq(min(ymz),max(ymz),median(ymz)*5*10^-6))
    ic=as.numeric(cut(ieic$y,unique(c(0,minHeight1,minHeight2,Inf))))
    l=which(ieic$y>=minHeight1)
    points(ieic$mz[l],ieic$rt[l],pch=16,col=ic[l],cex=c(.2,1,.8)[ic[l]])
    #      abline(v=range(ieic$mz,na.rm=TRUE))
  }
  lusid=unique(ieic$Sid[ieic$y>=minHeight2])
  ### loop over each sample
  allre=NULL
  for(ii in lusid){
    l2=which(ieic$y>minHeight1)
    llx<-.MGinquickSplit(ieic[l2,],ncons=parDeco$ncons,dppm = parDeco$ppm2,dmz=parDeco$dmz,minNoise=minHeight1)[[1]]
    if(length(llx)==0) next
    ldrmz=sapply(llx,function(x) diff(range(ieic[l2[x],"mz"])))
    l=which(ldrmz<=max(quantile(ldrmz,.75),mindmz))
    mllx=lapply(GRMeta:::.GRmergellx(llx[l]),function(x) l2[unique(x)])
    re=do.call("rbind",lapply(mllx,function(x) c(range(ieic$mz[x]),range(ieic$rt[x]),max(ieic$y[x]),ieic$mz[x[which.max(ieic$y[x])]])))
    #
    if(doPlot)  apply(re,1,function(x) rect(x[1],x[3],x[2],x[4]))
    allre=rbind(allre,re)
  }
  if(is.null(allre)) return(NULL)
  isover=.GRisover2(allre[,1],allre[,2],allre[,3]-addrt,allre[,4]+addrt,retOne = T,thr1 =mindmz,thr2=10^-4)
  newwin=do.call("rbind",lapply(isover,function(x){
    y=c(mz=range(allre[x,1:2]),rt=range(allre[x,3:4]),mzmed=matrixStats:::weightedMedian(allre[x,6],allre[x,5]))
    c(y,range(y[c(5,5,1,1,2,2)]*(1+0.5*c(-parDeco$ppm,parDeco$ppm,-1,1,-1,1)*10^-6),y[1:2]),int=max(allre[x,5]))
  }))
  newwin=newwin[which(newwin[,'int']>=minHeight2),,drop=F]
  if(nrow(newwin)==0) return(NULL)
  if(nrow(newwin)>1){
    isover=.GRisover2(newwin[,6],newwin[,7],newwin[,3]-addrt,newwin[,4]+addrt,retOne = T,thr1 =mindmz,thr2=10^-4)
    newwin=do.call("rbind",lapply(isover,function(x){
      c(range(newwin[x,1:2]),range(newwin[x,3:4]),
        mzmed=matrixStats:::weightedMedian(newwin[x,5],newwin[x,"int"]),range(newwin[x,6:7]),int=max(newwin[x,"int"]))
    }))
  }
  leics=list()
  for(i in 1:nrow(newwin)){
    x=newwin[i,]
    l=which(eic$rt>=(x[3]-addrt) & eic$rt<=(x[4]+addrt))
    neic=eic[l,]
    l=which(neic$mz<(x[6]) | neic$mz>x[7])
    if(length(l)>0){
      neic$mz[l]=neic$y[l]=NA
      ### sids without data
      sid2rm=names(which(tapply(is.na(neic$mz),neic$Sid,all)))
      if(length(sid2rm)>0) neic=neic[!neic$Sid%in%sid2rm,]
      ### sids with duplicated scans
      sid2dups=names(which(tapply(neic$scan,neic$Sid,function(x) max(table(x)))>1))
      if(length(sid2dups)>0) for(j in sid2dups){
        lnna=neic$scan[which(!is.na(neic$mz) & neic$Sid==j)]
        lna=neic$scan[which(is.na(neic$mz) & neic$Sid==j)]
        lsc2rm=intersect(lna,lnna)
        if(length(lsc2rm)>0) neic=neic[which(!(is.na(neic$mz) & neic$Sid==j & neic$scan%in%lsc2rm)),]
      }
    }  
    # ldups=which(duplicated(neic[,c("Sid","scan","id")]))
    # if(length(ldups)) stop('dupps')
    leics[[i]]=neic
  }
  # leics[[ix]]=NULL
  
  if(doPlot){
    apply(newwin,1,function(x) rect(x[6],x[3]-1*addrt,x[7],x[4]+1*addrt,border = 4,lwd=2))
    par(mfrow=c(1,1),mar=c(5.1 ,4.1 ,4.1, 2.1))
  }
  
  return(leics)
}
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 


######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
infctMRoi.get2<-function (fileres, what = "ROImat"){
  fileres = fileres[file.exists(fileres)]
  if (length(fileres) == 0) 
    stop("No file found")
  if (is.null(names(fileres))) 
    names(fileres) = paste0("S", 1:length(fileres))
  arois = list()
  for (isid in names(fileres)) {
    cat(".")
    load(fileres[isid])
    if(!exists(what)) next
    eval(parse(text=paste0("y=",what)))
    y$Sid = isid
    arois[[isid]] = y
    rm(list=what)
  }
  arois = do.call("rbind", arois)
  if (any(grepl("mz", names(arois)))) 
    arois = arois[order(rowMeans(arois[, grep("mz", names(arois))])), ]
  rownames(arois) = NULL
  return(arois)
}
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
