#' @name extractPI
#' @title Extract ROI around spectra
#' 
#' some stufff
#'
#' @param obj metaboGoS obj
#' @param parDeco deconvolution parameter list
#' @param minPIwindow window size in min. around the spectra
#' @param PIweight weighting
#' @param nSlaves number of slaves for parallel
#' @import foreach
#' @import doParallel
#' @return List of stuf
#' @export
extractPI<-function(obj,
                    parDeco=list(minPurity=c(0.2,0.4),
                         ppm=5,dmz=0.001,
                         drt=1,rtlim=c(-Inf,-Inf,Inf,Inf),
                         minHeightMS1=100000,minNoiseMS1=10000,minZero=200),
                    minPIwindow=NA,PIweight=2,nSlaves=1){
  
  ## minPIwindow=0.05;PIweight=2;nSlaves=10
  if(!"ddaSet"%in%class(obj)) stop("Error: not a ddaSet object")

  strt=as.integer(as.POSIXct( Sys.time() ))

  ms2inf=obj$MS2Infos
  xr=xcmsRaw(obj$File$File,includeMSn = FALSE)
  cat("Extract potential PIs from '",unclass(xr@filepath),"' comprising ",nrow(ms2inf)," MS/MS and ",length(xr@scantime)," MS1\n",sep="")
  
  sc2rt=round(xr@scantime/60,5)
  drt=quantile(diff(xr@scantime/60),.1)
  if(is.na(minPIwindow)) minPIwindow=quantile(diff(xr@scantime/60),9)+drt ## set it to near largest delta rt + drt
  
  ### Get potential precursor within minPIwindow of the spectra
  amzrt=do.call("rbind",lapply(1:nrow(ms2inf),function(ipmz) {
    scrange=range(which(abs(sc2rt-sc2rt[ms2inf[ipmz,"ScMS1"]])<=minPIwindow))
    re=GRMeta:::.GRrawMat(xr,scanrange=scrange,mzrange=ms2inf[ipmz,]$PrecMZ+c(-1,1)*(parDeco$winMS2+parDeco$dmz))
    if(nrow(re)==0) return(NULL)# return(c(PInt=NA,MZ=NA,SBR=NA,CodaWide=NA))
    we=.MGgetMZweight(re[,"mz"]-ms2inf[ipmz,]$PrecMZ,PIweight)
    re[,"y"]=re[,"y"]*we
    pint=re[,"y"]/tapply(re[,"y"],re[,"scan"],sum)[as.character(re[,"scan"])]
    if(max(pint)<parDeco$minPurity[1]) return(NULL)# return(c(PInt=NA,MZ=NA,SBR=NA,CodaWide=NA))
    re=cbind(i=ipmz,pint=pint,re[,c("mz","rt","y"),drop=F])
    re[pint>=parDeco$minPurity[1],,drop=F]}))
  cat(" -- ",nrow(amzrt)," ions found: ",length(unique(amzrt[,1]))," spectra out of ",nrow(ms2inf),"\n Purity summary:\n",sep="")
  print(summary(amzrt[,2]))
  
  
  llppm<-GRMeta:::.GRsplistMZ(amzrt[,"mz"],dppm = parDeco$ppm,dmz=parDeco$dmz,typ="min") ## min distance b/w dppm/dmz
  ROImat=do.call("rbind",lapply(llppm,function(x) c(median(amzrt[x,"mz"]),range(amzrt[x,"mz"]),round(range(amzrt[x,"rt"]),4))))
  colnames(ROImat)=c("mz","mzmin","mzmax","rtmin","rtmax")
  ROImat=ROImat[order(ROImat[,1]),,drop=F]
  cat(" -- reduced to ",nrow(ROImat)," ROIs [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")
  
  
  ### merge wide m/z/RT windows
  ROImat1=lapply(1:nrow(ROImat),function(ix){
    xroi=ROImat[ix,]
    lmz=range(c(xroi[c(1,1,2,3)]*(1+1*c(-2,2,-1,1)*parDeco$ppm*10^-6),xroi[2:3]+c(-1,1)*parDeco$dmz))+c(-1.1,1.1)*parDeco$dmz ## very large window!!
    lrt=xroi[4:5]+c(-1,1)*parDeco$drt
    xr=rawEIC(xr, mzrange =lmz)#, rtrange = lrt*60)
    l=which(xr$intensity>=parDeco$minZero)
    ll=GRMeta:::.GRsplist(sc2rt[l],l,d=parDeco$drt*2+2.1*drt)
    do.call("rbind",lapply(ll,function(x){ 
      newx=sc2rt[range(xr$scan[x])]+c(-1,1)*(parDeco$drt)
      #x2=x[xr$scan[x]>=rangeSc[2] & xr$scan[x]<=rangeSc[3]]
      #if(length(x2)<1) return(c(xroi[1],lmz[1],lmz[2],newx,0,0,0))
      aprt=xr$scan[x[which.max(xr$intensity[x])]]
      return(c(xroi[1],lmz[1],lmz[2],newx,diff(newx),max(xr$intensity[x]),sc2rt[aprt]))
    }))
  })
  ROImat1=do.call("rbind",ROImat1)
  colnames(ROImat1)=c(colnames(ROImat),"drt","intensity","aprt")
  
  l2k=which(ROImat1[,"intensity"]>=parDeco$minHeightMS1 & ROImat1[,"drt"]>=minPIwindow & ROImat1[,"rtmin"]<=parDeco$rtlim[4] & ROImat1[,"rtmax"] >= parDeco$rtlim[1])
  ROImat1=ROImat1[l2k,]
  ROImat1[ROImat1[,"rtmin"]<=parDeco$rtlim[1],"rtmin"]=parDeco$rtlim[1]
  ROImat1[ROImat1[,"rtmax"]>=parDeco$rtlim[4],"rtmax"]=parDeco$rtlim[4]
  
  ### Combine dupl -> remerge retention time -> could be improved!?!
  llover=GRMeta:::.GRisover(ROImat1[,"mzmin"],ROImat1[,"mzmax"],retOne = T)
  ROImat2=do.call("rbind",lapply(llover,function(x){
    if(length(x)==0) return(ROImat1[x,,drop=F])
    c(mean(ROImat1[x,"mz"]),min(ROImat1[x,"mzmin"]),max(ROImat1[x,"mzmax"]),min(ROImat1[x,"rtmin"]),max(ROImat1[x,"rtmax"]),0,
      max(ROImat1[x,"intensity"]),ROImat1[x[which.max(ROImat1[x,"intensity"])],"intensity"])
  }))
  colnames(ROImat2)=colnames(ROImat1)
  ROImat2[,"drt"]=ROImat2[,"rtmax"]-ROImat2[,"rtmin"]
  l2k=which(ROImat2[,"intensity"]>=parDeco$minHeightMS1 & ROImat2[,"drt"]>=minPIwindow & ROImat2[,"rtmin"]<=parDeco$rtlim[4] & ROImat2[,"rtmax"]>=parDeco$rtlim[1])
  ROImat2=ROImat2[l2k,]
  
  cat(" -- found: ",nrow(ROImat2)," potential ROIs [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")
  
  ####################### Refine frames
  llre=list()
  ll=1:nrow(ROImat2)
  
  if(nSlaves>1)   nSlaves=max(1, min(nSlaves,detectCores()-1))
  if(nSlaves>1){
    clProc<-makeCluster(nSlaves)
    registerDoParallel(clProc)
    cat(" -- registering ",nSlaves," clusters\n",sep="")
  }
  
   ### Parallele bit
  if(nSlaves>1) 
    #llre=foreach(i = ll, .export = fct2exp,.packages = c("igraph","xcms","GRMeta"), .verbose =F)  %dopar%{
      llre=foreach(i = ll,.packages = c("metaboGoS"), .verbose =F)  %dopar%{
        lmz=range(ROImat2[i,c("mzmin","mzmax")])
      lrt=range(ROImat2[i,c("rtmin","rtmax")])*60
      eic=GRMeta:::.GRrawMat(xr,mzrange = lmz, rtrange = lrt,padsc =F)
      if(nrow(eic)<2) return(list())
      if(max(eic[,"y"])<parDeco$minHeightMS1) return(list())
      re=.MGrefineROIs(eic,parDeco,minRTwin = minPIwindow,drt=drt)
      if(is.null(re)) return(list())
      cbind(id=i,re)
    }
  ## Serial bit
  lperc=round(seq(1,nrow(ROImat2),length=12)[2:11])
  if(nSlaves<=1) for(i in ll){
    if(i %in% lperc) cat(i," ")
    lmz=range(ROImat2[i,c("mzmin","mzmax")])
    lrt=range(ROImat2[i,c("rtmin","rtmax")])*60
    eic=GRMeta:::.GRrawMat(xr,mzrange = lmz, rtrange = lrt,padsc =F)
    if(nrow(eic)<2) next
    if(max(eic[,"y"])<parDeco$minHeightMS1) next
    re=.MGrefineROIs(eic,parDeco,drt=drt,minRTwin = minPIwindow,drt=drt)
    if(is.null(re)) next
    llre=c(llre,list(cbind(id=i,re)))
  }
  
  if(nSlaves>1) stopCluster(clProc)
  
  ## combine results
  ROImat3=do.call("rbind",llre[sapply(llre,length)>0])
  colnames(ROImat3)=c("id","rtmin","rtmax","mz50","mz10","mz90","intensity","mzap","mzmin","mzmax")
  l2k=which(ROImat3[,"intensity"]>=parDeco$minHeightMS1 & ROImat3[,"rtmin"]<=parDeco$rtlim[4] & ROImat3[,"rtmax"]>=parDeco$rtlim[1])
  ROImat3=ROImat3[l2k,,drop=F]
  ROImat3=ROImat3[order(ROImat3[,'mz50'],ROImat3[,"rtmin"]),,drop=F]
  cat(" -- number of ROIs after refinement: ",nrow(ROImat3)," [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")
  
  ### Associate MS/MS to ROI
  llsplit=split(1:nrow(amzrt), ceiling(seq_along(1:nrow(amzrt))/200))
  Roi2sp=do.call("rbind",lapply(llsplit,function(x){
    ddmz=outer(amzrt[x,"mz"],ROImat3[,"mzmin"],"-")>=(-parDeco$dmz/2) & outer(amzrt[x,"mz"],ROImat3[,"mzmax"],"-")<=parDeco$dmz/2
    ddrt=outer(amzrt[x,"rt"],ROImat3[,"rtmin"],"-")>=0 & outer(amzrt[x,"rt"],ROImat3[,"rtmax"],"-")<=0
    re=which(ddrt & ddmz,arr=T)
    re[,1]=x[re[,1]]
    re}))
  Roi2sp=cbind(amzrt[Roi2sp[,1],],Roi2sp)
  Roi2sp=cbind(Roi2sp,rtms2=obj$MS2Infos[Roi2sp[,"i"],]$RT)
  lso=order(Roi2sp[,"i"],Roi2sp[,"col"],abs(Roi2sp[,"rt"]-Roi2sp[,"rtms2"]))
  Roi2sp=Roi2sp[lso,]
  ## keep the maximum Pint around the original RT
  llsp=unlist(tapply(1:nrow(Roi2sp),Roi2sp[,"i"],function(x) tapply(x,Roi2sp[x,"col"],function(y) y[1:min(3,length(y))])),recursive = F)
  llsp=sapply(llsp,function(x) x[which.max(Roi2sp[x,"pint"])])
  # llsp=unlist(tapply(1:nrow(Roi2sp),Roi2sp[,1],function(x) tapply(x,Roi2sp[x,6],function(y) y)),recursive = F)
  Roi2sp=Roi2sp[llsp,]
  Roi2sp=Roi2sp[Roi2sp[,"col"]%in%names(which(tapply(Roi2sp[,"pint"],Roi2sp[,"col"],max)>=parDeco$minPurity[2])),]
  lurois=sort(unique(Roi2sp[,"col"]))
  lnurois=1:length(lurois)
  Roi2sp[,"col"]=lnurois[match(Roi2sp[,"col"],lurois)]
  ROImat4=ROImat3[lurois,-1]
  rownames(ROImat4)=sprintf("R%.4f@%.1f-%.1f",ROImat4[,"mz50"],ROImat4[,"rtmin"],ROImat4[,"rtmax"])
  
  SpId2ROI=data.frame(SpId=ms2inf$SpId[Roi2sp[,"i"]],
                     ROI=rownames(ROImat4)[Roi2sp[,"col"]],
                     PInt=round(Roi2sp[,"pint"],4),
                     yPrec=round(Roi2sp[,"y"],6),
                     mzPrec=round(Roi2sp[,"mz"],6),
                     rtPrec=round(Roi2sp[,"rt"],6))

  cat(" --> Final number of ROIs: ",nrow(ROImat4)," assoc. to ",length(unique(SpId2ROI$SpId)),"/",nrow(ms2inf),
      " MS/MS [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")

    
  invisible(list(ROI=ROImat4,Sp2ROI=SpId2ROI))
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## Refine ROIs
#' 
#' Refining ROIs
#' 
#' @keywords internal
#' 
#' @export
.MGrefineROIs<-function(eic,parDeco,minRTwin,drt){
 
  
  ll2use=GRMeta:::.GRsplitBoth2(eic[,"mz"],eic[,"rt"],dppm=parDeco$ppm*2,dmz=parDeco$dmz*2,drt=parDeco$drt+7*drt) ## allow +/- 3 pseudo scan around window
  ## only keep segments 
  todo=which(sapply(ll2use,function(x){
    l2chk=x[which(eic[x,"rt"]>=parDeco$rtlim[2] & eic[x,"rt"]<=parDeco$rtlim[3])]
    diffrt=ifelse(length(l2chk)>1,diff(range(eic[l2chk,"rt"])),0)
    any(eic[l2chk,"y"]>parDeco$minHeightMS1) & diffrt>=minRTwin}))
  
  if(length(todo)==0) return(NULL)
  ll2use=ll2use[todo]
  winmat=NULL
  ## loop over segments
  if(length(ll2use)) for(idx in 1:length(ll2use)){
    xeic=eic[ll2use[[idx]],]
    rmz=range(xeic[,"mz"])
    mindmz=min(c(mean(rmz)*parDeco$ppm*10^-6,parDeco$dmz))/2
    
    re=.MGgetConsecutive(xeic,ncons=2,ppm=parDeco$ppm,dmz=parDeco$dmz,mindmz,minRTwin=minRTwin,minNoise=parDeco$minNoiseMS1,minHeight=0)
      
    if(is.null(re[[1]])) next
    l=which(re$st[,6]>=parDeco$minHeightMS1)
    if(length(l)==0) next
    winmat=rbind(winmat,.MGgetWindow(re$ll[l],xeic,ppm=parDeco$ppm,dmz=parDeco$dmz,minRTwin=minRTwin,maxRTwin=parDeco$drt/2,maxIter = 6))
#    colnames(winmat)
  }
  return(winmat)
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### looking for consecutive scans
#' 
#' Getting consecutive m/z within ppm/dmz ROIs
#' 
#' @keywords internal
#' 
#' @export
.MGgetConsecutive<-function(xeic,ncons=3,ppm=5,dmz=0.001,mindmz=dmz/2,minRTwin=0.01,minNoise=10000,minHeight=minNoise*2){
  
  llx<-.MGinquickSplit(xeic,ncons = ncons,dppm = ppm,dmz=dmz,minNoise = minNoise)
  if(length(llx[[1]])==0) return(list(NULL,NULL))
  
  ldrmz=sapply(llx[[1]],function(x) diff(range(xeic[x,"mz"])))
  l=which(ldrmz<=max(quantile(ldrmz,.75),mindmz))
  
  # cols=(1:nrow(xeic)%in%unique(unlist(llx[[1]])))*1+(1:nrow(xeic)%in%unique(unlist(llx[[1]][l])))*1+1
  # if(doPlot) plot(xeic,cex=.5,pch=16,col=cols,main=imain)
  
  mllx=lapply(GRMeta:::.GRmergellx(llx[[1]][l]),unique)
  
  llx2=llx[[1]][-l]
  add2close=do.call("cbind",lapply(llx2,function(y) sapply(mllx,function(x) any(y%in%x))))
  toloop=ifelse(length(llx2)>0,any(colSums(add2close)==1),FALSE)
  while(toloop){
    l= which(colSums(add2close)==1)
    for(i in l) mllx[[which(add2close[,i])]]=c(mllx[[which(add2close[,i])]],llx2[[i]])
    llx2=llx2[-l]
    if(length(llx2)==0) break
    add2close=do.call("cbind",lapply(llx2,function(y) sapply(mllx,function(x) any(y%in%x))))
    toloop=ifelse(length(llx2)>0,any(colSums(add2close)==1),FALSE)
  }
  mllx=lapply(mllx,function(x) sort(unique(x)))
  # if(doPlot) for(i in 1:length(mllx)) lines(xeic[mllx[[i]],1:2])
  
  
  if(length(llx2)>0){ 
    llx2=GRMeta:::.GRmergellx(llx2)
    add2close=do.call("cbind",lapply(llx2,function(y) sapply(mllx,function(x) any(y%in%x))))
    for(ik in 1:2){
      toloop=ifelse(length(llx2)>0,any(colSums(add2close)==ik),FALSE)
      while(toloop){
        l= which(colSums(add2close)==ik)
        for(i in l) for(k in which(add2close[,i])) mllx[[k]]=c(mllx[[k]],llx2[[i]])
        llx2=llx2[-l]
        mllx=GRMeta:::.GRmergellx(mllx)
        if(length(llx2)==0) break
        add2close=do.call("cbind",lapply(llx2,function(y) sapply(mllx,function(x) any(y%in%x))))
        toloop=ifelse(length(llx2)>0,any(colSums(add2close)==ik),FALSE)
      }
    }
  }
  if(length(llx2)>0) llx2=llx2[sapply(llx2,function(x) max(xeic[x,3]))>minHeight]
  
  #if(doPlot) if(length(llx2)>0)  for(i in 1:length(llx2))  lines(xeic[llx2[[i]],1:2],lwd=2,col=2)
  if(length(llx2)>0) mllx=c(mllx,llx2)
  
  stllx=do.call("rbind",lapply(mllx,function(x) c(range(xeic[x,"rt"]),quantile(xeic[x,"mz"],c(.5,.1,.9)),max(xeic[x,"y"]),
                                                  xeic[x[which.max(xeic[x,"y"])],"mz"],range(xeic[x,"mz"]))))
  lso=order(stllx[,3])
  stllx=stllx[lso,,drop=F]
  mllx=mllx[lso]
  
  for(ik in 1:3){
    winsize=(ik+.1)*minRTwin
    if(nrow(stllx)>1)  hasovermz<-GRMeta:::.GRisover2(stllx[,4],stllx[,5],stllx[,1],stllx[,2],retOne = T,thr1=mindmz,thr2=winsize) else hasovermz=list(1)
    while(any(sapply(hasovermz,length)>1)){
      mllx<-lapply(hasovermz,function(x) sort(unique(unlist(mllx[x]))))
      stllx=do.call("rbind",lapply(mllx,function(x) c(range(xeic[x,"rt"]),quantile(xeic[x,"mz"],c(.5,.1,.9)),max(xeic[x,"y"]),
                                                      xeic[x[which.max(xeic[x,"y"])],"mz"],range(xeic[x,"mz"]))))
      if(nrow(stllx)>1) hasovermz<-GRMeta:::.GRisover2(stllx[,4],stllx[,5],stllx[,1],stllx[,2],retOne = T,thr1=mindmz,thr2=winsize) else hasovermz=list(1)
    }
  }
  
  #### fill in gaps
  for(i in 1:nrow(stllx)){
    x=stllx[i,]
    l=which(xeic[,"rt"]>=x[1] & xeic[,"rt"]<=x[2] & 
              ((xeic[,"mz"]>=x[4] & xeic[,"mz"]<=x[5]) | abs(1-xeic[,"mz"]/x["mz"])<ppm*10^-6  | abs(1-xeic[,"mz"]/x[3])<ppm*10^-6))
    mllx[[i]]=sort(unique(c(mllx[[i]],l)))
  }
  stllx=do.call("rbind",lapply(mllx,function(x) c(range(xeic[x,"rt"]),quantile(xeic[x,"mz"],c(.5,.1,.9)),max(xeic[x,"y"]),
                                                  xeic[x[which.max(xeic[x,"y"])],"mz"],range(xeic[x,"mz"]))))
  
  lso=order(stllx[,3])
  stllx=stllx[lso,,drop=F]
  mllx=mllx[lso]
  
  return(list(ll=mllx,st=stllx))
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### expand/reduces windows around initial guesses
#' 
#' Merge consecutive set of scans
#' 
#' @param mllx initial list of cons scans
#' @param xeic scan/mz/rt/y matrix
#' @param ppm tolerance in ppm
#' @param dmz tolerance in Da
#' @param minRTwin window stuff
#' @param maxRTwin window stuff
#' @param maxIter max number of iteration, usually done within 1-3
#' 
#' @keywords internal
#' @export
.MGgetWindow<-function(mllx,xeic,ppm,dmz,minRTwin,maxRTwin,maxIter=6){

  
  ## Comp stats
  stllx=stllx0=do.call("rbind",lapply(mllx,function(x) c(range(xeic[x,"rt"]),quantile(xeic[x,"mz"],c(.5,.1,.9)),max(xeic[x,"y"]),
                                                         xeic[x[which.max(xeic[x,"y"])],"mz"],range(xeic[x,"mz"]))))
  
  lmiss0=lmiss=1:nrow(xeic)
  lmiss=lmiss[!lmiss%in%unlist(mllx)]
  
  if(length(lmiss)==0 & nrow(stllx)==1) return(stllx)
  iround=1
  winover=maxdrt=minRTwin*3.1
  rmz=range(xeic[,"mz"])
  mindmz=min(c(mean(rmz)*ppm*10^-6,dmz))/2
  
  doLoop=TRUE
  while(doLoop){
    
    lmiss=lmiss[!lmiss%in%unlist(mllx)]
    
    mdppm=abs(10^6*(1-outer(xeic[lmiss,"mz"],stllx[,7],"/")))
    mdmz=abs(outer(xeic[lmiss,"mz"],stllx[,7],"-"))
    mdppm50=abs(10^6*(1-outer(xeic[lmiss,"mz"],stllx[,3],"/")))
    mdmz50=abs(outer(xeic[lmiss,"mz"],stllx[,3],"-"))
    mdmzwi=(mdppm<=ppm | mdmz<=dmz | mdppm50<=ppm | mdmz50<=dmz)
    mdmzin=outer(xeic[lmiss,"mz"],stllx[,4],"-")>=0 & outer(xeic[lmiss,"mz"],stllx[,5],"-")<=0 ## in the 10/90
    hasmz=(mdmzwi | mdmzin)
    
    mdrtle=outer(xeic[lmiss,"rt"],stllx[,1],"-") ## pos inside
    mdrtri=outer(xeic[lmiss,"rt"],stllx[,2],"-") ## neg inside
    
    ## check if anything to add
    toadd=which(hasmz & mdrtle>=(-maxdrt)  & mdrtri<=maxdrt,arr=T)
    if(length(toadd)) if(nrow(toadd)>0){
      toadd[,1]=lmiss[toadd[,1]]
      toadd=toadd[!duplicated(toadd[,1]),,drop=F]
      if(nrow(toadd)>0){
        for(k in unique(toadd[,2])) mllx[[k]]=c(mllx[[k]],toadd[toadd[,2]==k,1])
        stllx=do.call("rbind",lapply(mllx,function(x) c(range(xeic[x,"rt"]),quantile(xeic[x,"mz"],c(.5,.1,.9)),max(xeic[x,"y"]),
                                                        xeic[x[which.max(xeic[x,"y"])],"mz"],range(xeic[x,"mz"]))))
      }}
    
    ## merge windows -> greater window for each round
    if(nrow(stllx)>1)  
      hasovermz<-GRMeta:::.GRisover2(stllx[,4],stllx[,5],stllx[,1],stllx[,2],retOne = T,thr1=mindmz,thr2=winover) else    hasovermz=list(1)
    while(any(sapply(hasovermz,length)>1)){
      mllx<-lapply(hasovermz,function(x) sort(unique(unlist(mllx[x]))))
      stllx=do.call("rbind",lapply(mllx,function(x) c(range(xeic[x,"rt"]),quantile(xeic[x,"mz"],c(.5,.1,.9)),max(xeic[x,"y"]),
                                                      xeic[x[which.max(xeic[x,"y"])],"mz"],range(xeic[x,"mz"]))))
      if(nrow(stllx)>1) hasovermz<-GRMeta:::.GRisover2(stllx[,4],stllx[,5],stllx[,1],stllx[,2],retOne = T,thr1=mindmz,thr2=winover) else hasovermz=list(1)
    }
    iround=iround+1
    doLoop=(nrow(stllx0)!=nrow(stllx0) | all(lmiss0%in%lmiss)) & iround<=maxIter
    ### first round->small rt expansion!
    if(iround==2) doLoop=TRUE
    maxdrt=maxRTwin+minRTwin
    winover=minRTwin*5.1
    ###
    lmiss0=lmiss
    stllx0=stllx
  }
  
  return(stllx)
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## same of for the in getROIsToF2.R
#' 
#' Quick split/merge consecutive scan within ppm or dmz
#' 
#' @param dat scan/mz/rt/y matrix
#' @param ncons number of cons scans
#' @param dppm tolerance in ppm
#' @param dmz tolerance in Da
#' @param minNoise only use points greter than minNoise
#' @keywords internal
#' 
#' @export
.MGinquickSplit<-function(dat,ncons=4,dppm=5,dmz=0.005,minNoise=1000){
  # dat=xeic;ncons=parDeco$ncons;dppm=parDeco$ppm;minNoise=parDeco$minNoise;dmz=parDeco$dmz
  
  vscan=dat[,"scan"]
  vint=dat[,"y"]
  vmz=dat[,"mz"]
  
  
  lsup=which(vint>=minNoise)
  luscn=min(vscan[lsup]):max(vscan[lsup])
  if((diff(range(luscn))+1)<ncons) return(NULL)
  lst=luscn[1:(length(luscn)-ncons+1)]
  lend=luscn[(ncons):length(luscn)]
  
  llx=unlist(lapply(1:length(lst),function(i){
    lx=lsup[which(vscan[lsup]%in%c(lst[i]:lend[i]))]
    GRMeta:::.GRsplistMZ(vmz[lx],iv = lx,dppm = dppm*1.001,dmz=dmz) ##
  }),recursive = F)
  
  # cat("0: ",Sys.time()-strtx,"\n",sep="");strtx=Sys.time()
  # print(length(llx))
  if(length(llx)==0) return(NULL)
  llxnsc=t(sapply(llx,function(x){y=unique(vscan[x]);c(length(y),range(y),length(x))}))
  l2k=which(llxnsc[,1]>=ncons & llxnsc[,1]>(llxnsc[,3]-llxnsc[,2]))
  llx=llx[l2k]
  llxnsc=llxnsc[l2k,,drop=F]
  l2chk=which(llxnsc[,4]>llxnsc[,1])
  if(length(l2k)>0){
    dmz2=min(1.001*10^-6*dppm*median(vmz[lsup]),dmz)
    llx2=unlist(lapply(llx[l2chk],function(lx) GRMeta:::.GRsplist(vmz[lx],iv = lx,dmz2,ismass = F)),rec=F)
    llxnsc2=do.call("rbind",lapply(llx2,function(x){y=unique(vscan[x]);c(length(y),range(y),length(x))}))
    l2k=which(llxnsc2[,1]>=ncons & llxnsc2[,1]>(llxnsc2[,3]-llxnsc2[,2]))
    llx2=llx2[l2k]
    llxnsc2=llxnsc2[l2k,,drop=F]
    if(length(llx2)>0){llx=c(llx,llx2)[-l2chk];llxnsc=rbind(llxnsc,llxnsc2)[-l2chk,,drop=F]}
  }
  list(llx,llxnsc)
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

