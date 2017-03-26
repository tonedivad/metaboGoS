### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Weighting around the precursor window -> to be improved with larger windo etc...
#' 
#' @param eic vector of mz-mzPrec
#' @param parDeco weight
#' @param idx weight
#' @param doPlot weight
#' @keywords internal
#' 
#' @export
.MGrefFct0<-function(eic,parDeco,idx,doPlot=T){
  
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
    rtr=rtr+c(-1,1)*parDeco$psdrt*(parDeco$span+1)  ## extra padding
    xeic=eic[which(eic[,"rt"]>=rtr[1] & eic[,"rt"]<=rtr[2]),]
    l2excl=which(xeic[,"mz"]<mzr[1] | xeic[,"mz"]>mzr[2])
    if(length(l2excl)) xeic[l2excl,c("mz","y")]=NA
    xeic=xeic[order(xeic[,"scan"],!is.na(xeic[,"mz"])),]
    xeic=xeic[which(!(duplicated(xeic[,"scan"]) & is.na(xeic[,"mz"]))),]
    
    ire=.MGintegrateEIC(xeic,drt = parDeco$psdrt,span=parDeco$span,bw=parDeco$bw,minNoise=parDeco$minNoiseMS1,minHeight=parDeco$minHeightMS1,sbslr=2)
    if(is.null(ire[[1]])) next
    pks=ire[[1]]
    newx=ire[[2]]
    tmproi=rep(1,nrow(pks))
    if(nrow(pks)>1)
      for(i in 2:(nrow(pks))) tmproi[i]=tmproi[i-1]+ifelse((pks$pk.left[i]-pks$pk.right[i-1])>=3*parDeco$drt,1,0)
    
    llpks=tapply(1:nrow(pks),tmproi,function(x) pks[x,])
    

    ### reformat the ROI
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
      abline(h=parDeco$minHeightMS1)
      points(pks$pk.loc,pks$pk.int,pch=16,col=pks$pk.cl,cex=2)
      abline(v=pks$pk.left+parDeco$psdrt/2,lty=3)
      abline(v=pks$pk.right-parDeco$psdrt/2,lty=4)
    }
    
    reroi2=c(reroi2,list(pks))
  }
  reroi2=do.call("rbind",reroi2)
  return(reroi2)
}
