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
    
    ire=.MGintegrateEIC(xeic,drt = parDeco$psdrt,span=parDeco$span,bw=parDeco$bw,minNoise=parDeco$minNoiseMS1,minHeight=parDeco$minHeightMS1)
    if(is.null(ire[[1]])) next
    pks=ire[[1]]
    nrtr=range(c(pks$pk.left,pks$pk.right))+c(-1,1)*parDeco$drt+c(-1,1)*parDeco$psdrt*parDeco$span
    rtr[1]=max(rtr[1],nrtr[1])
    rtr[2]=min(rtr[2],nrtr[2])
    pks$rtmin=rtr[1]
    pks$rtmax=rtr[2]
    
    ### reformat the ROI
    l=which(xeic[,"rt"]>=rtr[1] & xeic[,"rt"]<=rtr[2] & !is.na(xeic[,"mz"]))
    ieic=xeic[l,]
    pks$mzmin=min(ieic[,"mz"]);pks$mzmax=max(ieic[,"mz"])
    pks$mz50=median(ieic[,"mz"]);pks$mz10=quantile(ieic[,"mz"],.1);pks$mz90=quantile(ieic[,"mz"],.9)
    roiap=which.max(ieic[,"y"])
    pks$intensity=ieic[roiap,"y"]
    pks$rt=ieic[roiap,"rt"]
    pks$mz=ieic[roiap,"mz"]
    pks$roi=idx
    pks$subroi=iiroi
    ###
    
    if(doPlot){
      newx=ire[[2]]
      
      plot(xeic[,"rt"],xeic[,"y"],pch=16,main=paste(idx,iiroi,round(.GRcoda(newx$y),2)),xlab="Retention time",ylab="Height",
           xlim=range(pretty(rtr)),ylim=range(c(0,xeic[,"y"],newx$y),na.rm=T),yaxt="n")
      axis(2,las=2)
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
