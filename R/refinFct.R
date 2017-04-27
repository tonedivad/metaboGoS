### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Weighting around the precursor window -> to be improved with larger windo etc...
#' TO DO: align properly the deconv/peak picking with the EIC integration!!!
#' 
#' @param eic vector of mz-mzPrec
#' @param eicbl vector of mz-mzPrec for the blank samples
#' @param parDeco weight
#' @param imain for the plot tile
#' @param doPlot weight
#' @keywords internal
#' 
#' @export
.MGrefFct0<-function(eic,eicbl=NULL,parDeco,imain=NULL,doPlot=T,lExcl=NULL){
  
  reroi=.MGrefineROIs(eic[which(!is.na(eic[,"y"])),],parDeco,minRTwin = parDeco$minRTwin,drt=parDeco$psdrt)
  if(is.null(reroi)) return(NULL)
  l2k=which(reroi[,"intensity"]>=parDeco$minHeightMS1 & (reroi[,"rtmax"]-reroi[,"rtmin"])>parDeco$minRTwin)
  if(length(l2k)==0) return(NULL)
  reroi=reroi[l2k,,drop=F]
  
  nspan=floor(parDeco$span/2)*2+1
  minNoise=parDeco$minNoiseMS1
  llmatpks=list()
  for(iiroi in 1:nrow(reroi)){
    
    #########
    iroi=reroi[iiroi,]
    mzr=range(iroi[c("mzmin","mzmax")])
    rtr0=rtr=range(iroi[c("rtmin","rtmax")])
    rtr=rtr+c(-3,3)*parDeco$psdrt*(parDeco$span+1)  ## extra padding
    xeic=eic[which(eic[,"rt"]>=rtr[1] & eic[,"rt"]<=rtr[2]),]
    l2excl=which(xeic[,"mz"]<mzr[1] | xeic[,"mz"]>mzr[2])
    if(length(l2excl)) xeic[l2excl,c("mz","y")]=NA
    xeic=xeic[order(xeic[,"scan"],!is.na(xeic[,"mz"])),]
    xeic=xeic[which(!(duplicated(xeic[,"scan"]) & is.na(xeic[,"mz"]))),]
    irmz=range(xeic[,"mz"],na.rm = T)
    
    ##### dual smooth/peak detec
    # drt = parDeco$psdrt ; span=parDeco$span;bw=parDeco$bw;minNoise=parDeco$minNoiseMS1;minHeight=parDeco$minHeightMS1;sbslr=2
    newx=.MGfill(xeic,parDeco$psdrt,nspan,parDeco$bw,minNoise)
    # y=newx[,"y2"]
    # n2pad=floor(parDeco$drt/parDeco$psdrt/2)
    # y=c(rep(y[1],n2pad),y,rep(rev(y)[1],n2pad))
    # bsl=GRMeta:::.GRbslrf(1:length(y),y,NoXP = NULL)
    # bsl$fit[bsl$fit<minNoise]=minNoise
    # bslscore <- (y - bsl$fit)/max(bsl$sigma, 10^-3)
    # bslscore[which(abs(bslscore) > 10)] = sign(bslscore[which(abs(bslscore) > 10)]) * 10
    # #newx$bslc=bslscore[nspan*5+1+(1:length(newx$y))]
    # newx=data.frame(cbind(newx,bsl=bsl$fit[n2pad+(1:nrow(newx))],bslc=bslscore[n2pad+(1:nrow(newx))]))
    # noise.local=(newx$y2/newx$bsl >=2 & newx$bslc>=1)*1;noise.local[newx$bslc< -1]=-1
    # 
    # bsl=newx$bsl;bslscore=newx$bslc;x=newx[,"rt"];y=newx[,"y2"];span=nspan;minNoise = minNoise*1.01;v2alim = 0.8;snr.thresh = parDeco$sbslr
    
    #### quick integration
    pks=.MGsimpleIntegr2(x=newx[,"rt"],y=newx[,"y2"],bsl = newx$bsl,bslscore = newx$bslc,snr.thresh = parDeco$sbslr,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8)
 #  pks=.MGsimpleIntegr(x=newx[,"rt"],y=newx[,"y2"],noise.local =newx$bslc,snr.thresh = 1,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8)
    if(nrow(pks)==0) next
    
    ##### check blank from and align to newx$x
    newx$ybl=rep(parDeco$minNoiseMS1,nrow(newx))
    if(!is.null(eicbl)){
      tmpbl=eicbl
      l2excl=which(tmpbl[,"mz"]<mzr[1] | tmpbl[,"mz"]>mzr[2])
      if(length(l2excl)) tmpbl[l2excl,c("mz","y")]=NA
      tmpbl=tmpbl[order(tmpbl[,"scan"],!is.na(tmpbl[,"mz"])),]
      tmpbl=tmpbl[which(!(duplicated(tmpbl[,"scan"]) & is.na(tmpbl[,"mz"]))),]
      tmpbl=tmpbl[tmpbl[,"blid"]%in%names(which(tapply(!is.na(tmpbl[,"y"]),tmpbl[,"blid"],sum)>10)),,drop=F]
      if(nrow(tmpbl)>1){
        tmpbl[which(is.na(tmpbl[,"y"]) | tmpbl[,"y"]<minNoise),"y"]=minNoise
        tmpbl<-try(do.call("cbind",tapply(1:nrow(tmpbl),tmpbl[,"blid"],function(xbl){
          ybl = .GRasysm(tmpbl[xbl,"y"], p = 0.5, lambda = 10^5)
          ybl[ybl <= minNoise | is.na(ybl)] = minNoise
          ybl=approx(tmpbl[xbl,"rt"],ybl,newx$rt,yleft = ybl[1],yright = rev(ybl)[1])$y})),TRUE)
        if(!"try-error"%in%class(tmpbl)) newx$ybl=apply(tmpbl,1,median)
      }
    }
    
    ## compute the SBR around the apex
    pks$pk.sbr=round(sapply(pks$pk.loc,function(i) max((newx$y2/newx$ybl)[which(abs(newx$rt-i)<=(2*parDeco$psdrt))])),3)
    pks$pk.snr=round(sapply(pks$pk.loc,function(i) max((newx$y2/newx$bsl)[which.min(abs(newx$rt-i))])),3)
    
    sbrlim=ifelse(any(lExcl>=irmz[1] & lExcl<=irmz[2]),0,parDeco$sbr)
    l2k=as.numeric(names(which(tapply(pks$pk.sbr,pks$pk.cl,max)>=sbrlim &
                                 tapply(pks$pk.snr,pks$pk.cl,max)>=parDeco$sbslr &
                                 tapply(pks$pk.int*2,pks$pk.cl,max)>parDeco$minHeightMS1)))
    pks=pks[pks$pk.cl%in%l2k,,drop=F]
    if(nrow(pks)==0) next
    pks$pk.cl=as.numeric(factor(pks$pk.cl))
    
    ## compute coda/roi stats
    pks$pk.roi=rep(1,nrow(pks))
    if(nrow(pks)>1)
      for(i in 2:(nrow(pks))) pks$pk.roi[i]=pks$pk.roi[i-1]+ifelse((pks$pk.left[i]-pks$pk.right[i-1])>=parDeco$drt,1,0)
    
    pks=cbind(pks,do.call("rbind",tapply(1:nrow(pks),pks$pk.roi,function(x){
      nrtr0=range(c(pks$pk.left[x],pks$pk.right[x]))+c(-.49,.49)*parDeco$drt
      nrtr0[1]=max(nrtr0[1],min(rtr0))
      nrtr0[2]=min(nrtr0[2],max(rtr0))
      l1=which(newx$rt>=nrtr0[1] & newx$rt<=nrtr0[2])
      nrtr=range(c(pks$pk.left[x],pks$pk.right[x]))
      l2=which(newx$rt>=nrtr[1] & newx$rt<=nrtr[2])
      data.frame(roi.id=sprintf("R%.4f@%.2f-%.2f",median(newx[l2,"mz"],na.rm=T),min(newx[l2,"rt"],na.rm=T),max(newx[l2,"rt"],na.rm=T)),
                 roi.coda=.GRcoda(newx$y2[l1]),roi.mz50=unname(median(newx[l2,"mz"],na.rm=T)),
                 roi.mz10=unname(quantile(newx[l2,"mz"],.1,na.rm=T)),roi.mz90=unname(quantile(newx[l2,"mz"],.9,na.rm=T)))
    }))[as.character(pks$pk.roi),,drop=FALSE])
    
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
                               pk.mzmin=round(min(tmp$mz[iap],na.rm = T),6),
                               pk.mzmax=round(max(tmp$mz[iap],na.rm = T),6),
                               pk.mzmed=suppressWarnings(round(matrixStats:::weightedMedian(tmp$mz,tmp$y,na.rm=T),6)),
                               pk.area.sm=round(unname(.GRgetArea(tmp$rt,tmp$y2)[1]),3),
                               pk.area=round(.GRgetArea(tmp$rt[!is.na(tmp$y)],tmp$y[!is.na(tmp$y)])[1],3),
                               pk.snr=tmp$y[iap]/tmp$bsl[iap])
    }
    if(length(matpks)==0) next
    matpks=do.call("rbind",matpks)
    
    matpks=cbind(matpks,pk.sbr=pks[matpks$id,c("pk.sbr")],pks[matpks$id,grep("^roi.",names(pks))])
    matpks$id=paste0("R",pks$pk.roi,"C",pks$pk.cl)[matpks$id]
    matpks$cl.id=tapply(1:nrow(matpks),matpks$id,function(x) sprintf("C%.5f@%.3f",median(matpks$pk.mz[x]),median(matpks$pk.rtap[x])))[matpks$id]
    llmatpks[[iiroi]]=matpks
    ###
    
    if(doPlot){
      #  pks=ire[[1]]
      # newx=ire[[2]]
      yl=pretty(range(sqrt(c(0,xeic[,"y"],newx$y)),na.rm=T))
     # yl=pretty(sqrt(range(pretty(yl^2))))
      xl=pretty(rtr0,n=11)
      par(mfrow=c(2,1))
      par(mar=c(1,5,2,0))
      plot(xeic[,"rt"],sqrt(xeic[,"y"]),pch=16,main=paste(imain,iiroi,round(.GRcoda(newx$y2),2)),xlab="",ylab="Height*10^6",
           xlim=range(xl),ylim=range(yl),axes=F)
      axis(2,at=yl,yl^2/1000000,las=2,pos=min(xl))
    #  axis(1,at=xl,pos=min(yl))
#      abline(v=c(unique(matpks$pk.rtmax),unique(matpks$pk.rtmin)),col="grey",lwd=2)
      lines(newx$rt,sqrt(newx$y2),col=2)
      lines(newx$rt,sqrt(newx$bsl),col=4,lwd=2,lty=1)
      lines(newx$rt,sqrt(newx$bsl*parDeco$sbslr),col=4,lwd=1,lty=4)
      lines(newx$rt,sqrt(newx$ybl*parDeco$sbr),col=3,lwd=1,lty=2)
      lines(newx$rt,sqrt(newx$ybl),col=3,lwd=2,lty=1)
      abline(h=sqrt(parDeco$minHeightMS1))
      points(matpks$pk.rtap,sqrt(matpks$pk.int.ap),pch=8,col=factor(matpks$cl.id),cex=2)
      abline(v=matpks$pk.rtmin+parDeco$psdrt/2,lty=3)
      abline(v=matpks$pk.rtmax-parDeco$psdrt/2,lty=4)
      
      
      ################
      par(mar=c(4,5,0,0))
      yl=pretty(range(c(range(matpks$pk.mzmed,matpks$pk.mz)*(1+c(-1,1)*2.1*10^-6),quantile(xeic[,"mz"],c(.1,.9),na.rm=T))))
      xl=pretty(rtr0,n=11)
      plot(xeic[,"rt"],xeic[,"mz"],pch=16,main=NULL,xlab="Retention time",ylab="",xlim=range(xl),ylim=range(yl),axes=F)
      axis(2,at=yl,las=2,pos=min(xl))
      nv=ceiling(diff(range(yl))/(median(yl)*10^-6))
      abline(h=(-nv:nv)*median(yl)*10^-6+median(yl),lty=3,col="grey")
      axis(1,at=xl,pos=min(yl))
#      abline(v=c(unique(matpks$pk.rtmax),unique(matpks$pk.rtmin)),col="grey",lwd=2)
      points(matpks$pk.rtap,matpks$pk.mz,pch=19,col=factor(matpks$cl.id),cex=2)
      abline(v=matpks$pk.rtmin+parDeco$psdrt/2,lty=3)
      abline(v=matpks$pk.rtmax-parDeco$psdrt/2,lty=4)
      par(mfrow=c(1,1),mar=c(5.1 ,4.1 ,4.1, 2.1))
      
      
      
    }
    
    #   reroi2=c(reroi2,list(pks))
  }
  allmatpks=do.call("rbind",llmatpks)
  return(allmatpks)
}
