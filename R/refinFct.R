### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Weighting around the precursor window -> to be improved with larger windo etc...
#' TO DO: align properly the deconv/peak picking with the EIC integration!!!
#' 
#' @param eic vector of mz-mzPrec
#' @param eicbl vector of mz-mzPrec for the blank samples
#' @param parDeco weight
#' @param idx weight
#' @param doPlot weight
#' @keywords internal
#' 
#' @export
.MGrefFct0<-function(eic,eicbl=NULL,parDeco,idx,doPlot=T){
  
  reroi=.MGrefineROIs(eic[which(!is.na(eic[,"y"])),],parDeco,minRTwin = parDeco$minRTwin,drt=parDeco$psdrt)
  if(is.null(reroi)) return(NULL)
  l2k=which(reroi[,"intensity"]>=parDeco$minHeightMS1 & (reroi[,"rtmax"]-reroi[,"rtmin"])>parDeco$minRTwin)
  if(length(l2k)==0) return(NULL)
  reroi=reroi[l2k,,drop=F]
  reroi2=list()
  for(iiroi in 1:nrow(reroi)){
    iroi=reroi[iiroi,]
    mzr=range(iroi[c("mzmin","mzmax")])
    rtr0=rtr=range(iroi[c("rtmin","rtmax")])
    rtr=rtr+c(-3,3)*parDeco$psdrt*(parDeco$span+1)  ## extra padding
    xeic=eic[which(eic[,"rt"]>=rtr[1] & eic[,"rt"]<=rtr[2]),]
    l2excl=which(xeic[,"mz"]<mzr[1] | xeic[,"mz"]>mzr[2])
    if(length(l2excl)) xeic[l2excl,c("mz","y")]=NA
    xeic=xeic[order(xeic[,"scan"],!is.na(xeic[,"mz"])),]
    xeic=xeic[which(!(duplicated(xeic[,"scan"]) & is.na(xeic[,"mz"]))),]
 
    ##### dual smooth/peak detec
    ire=.MGintegrateEIC(xeic,drt = parDeco$psdrt,span=parDeco$span,bw=parDeco$bw,minNoise=parDeco$minNoiseMS1,minHeight=parDeco$minHeightMS1,sbslr=2)
    if(is.null(ire[[1]])) next
    pks=ire[[1]]
    newx=ire[[2]]
    
    ##### check blank from and align to newx$x
    newx$ybl=rep(parDeco$minNoiseMS1,length(newx$x))
    if(!is.null(eicbl)){
      tmpbl=eicbl
      l2excl=which(tmpbl[,"mz"]<mzr[1] | tmpbl[,"mz"]>mzr[2])
      if(length(l2excl)) tmpbl[l2excl,c("mz","y")]=NA
      tmpbl=tmpbl[order(tmpbl[,"scan"],!is.na(tmpbl[,"mz"])),]
      tmpbl=tmpbl[which(!(duplicated(tmpbl[,"scan"]) & is.na(tmpbl[,"mz"]))),]
      tmpbl=tmpbl[tmpbl[,"blid"]%in%names(which(tapply(!is.na(tmpbl[,"y"]),tmpbl[,"blid"],sum)>10)),,drop=F]
      if(nrow(tmpbl)>1){
        tmpbl[which(is.na(tmpbl[,"y"]) | tmpbl[,"y"]<parDeco$minNoiseMS1),"y"]=parDeco$minNoiseMS1
        tmpbl<-try(do.call("cbind",tapply(1:nrow(tmpbl),tmpbl[,"blid"],function(xbl){
        ybl = .GRasysm(tmpbl[xbl,"y"], p = 0.5, lambda = 10^5)
        ybl[ybl <= parDeco$minNoiseMS1 | is.na(ybl)] = parDeco$minNoiseMS1
        ybl=approx(tmpbl[xbl,"rt"],ybl,newx$x,yleft = ybl[1],yright = rev(ybl)[1])$y})),TRUE)
      if(!"try-error"%in%class(tmpbl)) newx$ybl=apply(tmpbl,1,min)
      }
    }
    newx$filter=.GRfiltreScan((newx$y/newx$ybl)>2)
    ire[[2]]$ybl=newx$ybl
    
    ## only keep peaks above blank
    l2k=which(sapply(pks$pk.loc,function(i) any(newx$filter[which(abs(newx$x-i)<=(2*parDeco$psdrt))])))
    pks=pks[l2k,,drop=F]
    if(nrow(pks)==0) next
    ##### reformat rois : groups of peaks separated by at leat 3*drt
    tmproi=rep(1,nrow(pks))
    if(nrow(pks)>1)
      for(i in 2:(nrow(pks))) tmproi[i]=tmproi[i-1]+ifelse((pks$pk.left[i]-pks$pk.right[i-1])>=3*parDeco$drt,1,0)
    
    llpks=tapply(1:nrow(pks),tmproi,function(x) pks[x,])
    for(ipks in 1:length(llpks)){
      pks=llpks[[ipks]]
      rtr=rtr0
    nrtr=range(c(pks$pk.left,pks$pk.right))+c(-1.49,1.49)*parDeco$drt#+c(-1,1)*parDeco$psdrt*parDeco$span
    rtr[1]=max(rtr[1],nrtr[1])
    rtr[2]=min(rtr[2],nrtr[2])
    pks$rtmin=rtr[1]
    pks$rtmax=rtr[2]
    
    l=which(newx$x>=rtr[1] & newx$x<=rtr[2])
    icoda=.GRcoda(newx$y)
    l=which(xeic[,"rt"]>=rtr[1] & xeic[,"rt"]<=rtr[2] & !is.na(xeic[,"mz"]))
    ieic=xeic[l,]
    pks$mzmin=min(ieic[,"mz"]);pks$mzmax=max(ieic[,"mz"])
    pks$mz50=median(ieic[,"mz"]);pks$mz10=quantile(ieic[,"mz"],.1);pks$mz90=quantile(ieic[,"mz"],.9)
    roiap=which.max(ieic[,"y"])
    pks$intensity=ieic[roiap,"y"]
    pks$rt=ieic[roiap,"rt"]
    pks$mz=ieic[roiap,"mz"]
    pks$coda=icoda
    pks$roi=idx
    pks$subroi=iiroi
    pks$pk.cl=as.numeric(factor(pks$pk.cl))
    pks$newroi=sprintf("R%.4f@%.2f-%.2f",pks$mz50,pks$rtmin,pks$rtmax)
    llpks[[ipks]]=pks
    }
    pks=do.call("rbind",llpks)
    ###
    
    if(doPlot){
    #  pks=ire[[1]]
      newx=ire[[2]]
      
      plot(xeic[,"rt"],xeic[,"y"],pch=16,main=paste(idx,iiroi,round(.GRcoda(newx$y),2)),xlab="Retention time",ylab="Height",
           xlim=range(pretty(rtr0)),ylim=range(c(0,xeic[,"y"],newx$y),na.rm=T),yaxt="n")
      axis(2,las=2)
      abline(v=c(unique(pks$rtmax),unique(pks$rtmin)),col="grey",lwd=2)
      lines(newx,col=2)
      lines(newx$x,newx$bsl,col=4)
      lines(newx$x,newx$ybl,col=3)
      abline(h=parDeco$minHeightMS1)
      points(pks$pk.loc,pks$pk.int,pch=8,col=pks$pk.cl,cex=2)
      abline(v=pks$pk.left+parDeco$psdrt/2,lty=3)
      abline(v=pks$pk.right-parDeco$psdrt/2,lty=4)
    }
    
    reroi2=c(reroi2,list(pks))
  }
  reroi2=do.call("rbind",reroi2)
  return(reroi2)
}
