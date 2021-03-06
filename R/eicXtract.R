### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name eicXtract
#' 
#' Fetch ROIs/EICs across files 
#'
#' @param infile metaboGoS obj
#' @param eicmat eicmat
#' @param addppm add ppm 
#' @param refFct refining function
#' @param drt plotting at each step
#' @param nSlaves number of slaves for parallel
#' @param save numbl
#' @param outfile n
#' @param paddRT extra padding for RT
#' @import foreach
#' @import doParallel
#' @return metaboGoS obj with RoiInfos / PeakInfos /params
#' @export
eicXtract<-function(infile,eicmat,addppm=NA,refFct=metaboGoS:::.MGfill,drt,nSlaves=1,save=FALSE,outfile=NULL,paddRT=drt,...){
  
  cat("Processing ", infile, " - ",nrow(eicmat)," Roi/EICs to fecth\n",sep="")
  #drt=0.00446
  xr=xcmsRaw(infile)
  if(!is.na(addppm)) xr=.GRcorrRawPPM(xr,addppm)
  newrtx=(1:ceiling(max(xr@scantime/60/drt)))*drt
  sc2nrt=apply(abs(outer(newrtx,xr@scantime/60,"-")),2,which.min)
  sc2nrt=cbind(scan=1:length(xr@scantime),nscan=sc2nrt,nrt=round(newrtx[sc2nrt],5))
  #sc2nrt=round(newrtx[sc2nrt],6)
  
  allxeic=list()
  ll=eicmat$RoiId#[101:150]
  lperc=""
  if(length(ll)>20) lperc=ll[round(seq(1,length(ll),length=12)[2:11])]
  
  if(nSlaves>1)   nSlaves=max(1, min(nSlaves,detectCores()-1))
  if(nSlaves>1){
    clProc<-makeForkCluster(nSlaves)
    doParallel:::registerDoParallel(clProc)
    cat(" -- registering ",nSlaves," clusters\n",sep="")
  }
  
  ### Parallele bit
  if(nSlaves>1){
    #llre=foreach(i = ll, .export = fct2exp,.packages = c("igraph","xcms","GRMeta"), .verbose =F)  %dopar%{
    allxeic=foreach(iroi = ll,.packages = c("metaboGoS"), .verbose =F)  %dopar%{
      i=which(eicmat$RoiId==iroi)
      rtrange=(range(eicmat[i,c("rtmin","rtmax")])+c(-1,1)*paddRT)
      xeic=.GRrawMat(xr,mzrange=range(eicmat[i,c("mzmin","mzmax")]),rtrange = rtrange*60,padsc =TRUE) # extrat pad for the newscan
      if(all(is.na(xeic[,"y"]))) return(NULL)
      xeic=cbind(xeic,sc2nrt[xeic[,"scan"],2:3])
      if(!is.null(refFct)) df=data.frame(do.call("cbind",.MGfill(xeic,drt=drt,...))) else df=xeic
      return(list(iroi,df))
    }
    allxeic=allxeic[sapply(allxeic,length)==2]
    nallxeic=sapply(allxeic,function(x) x[[1]]) 
    allxeic=lapply(allxeic,function(x) x[[2]])
    names(allxeic)=nallxeic
  }
  ## Serial bit
  if(nSlaves<=1) for(iroi in ll){
    if(iroi %in% lperc) cat(iroi,"(",which(ll==iroi),") ",sep="")
    i=which(eicmat$RoiId==iroi)
    rtrange=(range(eicmat[i,c("rtmin","rtmax")])+c(-1,1)*paddRT)
    xeic=.GRrawMat(xr,mzrange=range(eicmat[i,c("mzmin","mzmax")]),rtrange = rtrange*60,padsc =TRUE) # extrat pad for the newscan
    if(all(is.na(xeic[,"y"]))) next
    xeic=cbind(xeic,sc2nrt[xeic[,"scan"],2:3])
    #if(any(table(xeic[,"scan"])>1)) print(c(sum(table(xeic[,"scan"])>1),i))
    if(!is.null(refFct)) allxeic[[iroi]]=data.frame(do.call("cbind",.MGfill(xeic,drt=drt,...))) else allxeic[[iroi]]=xeic
    
  }
  cat('---> Found ',length(allxeic),"/",nrow(eicmat),"\n",sep="")
  if(nSlaves>1) stopCluster(clProc)
  if(save & is.null(outfile)) save(file=paste0(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", infile),".rda"),allxeic,sc2nrt)
  if(!is.null(outfile)) save(file=outfile,allxeic,sc2nrt)
  invisible(allxeic)
  
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Merge ROIs/EICs from several files
#' 
#' @param files list of files
#' @param eics ROI/EICs names to recoup across files
#' 
#' @export
eicGather<-function(files,eics,conv2df=FALSE){
  
  
  files=files[file.exists(files)]
  if(length(files)==0) stop('No file found')
  if(is.null(names(files))){
    cat(" -- naming sids")
    names(files)=paste0("S",1:length(files))
  }
  
  ascnrt=list()
  aeics=vector("list",length(eics))
  names(aeics)=eics
  for(isid in names(files)){
    load(files[isid])
    allxeic=allxeic[names(allxeic)%in%eics]
    cat(isid," ",length(allxeic),"/",length(eics)," ",sep="")
    ascnrt[[isid]]=sc2nrt
    if(length(allxeic)==0) next
    for(i in names(allxeic)) aeics[[i]][[isid]]=allxeic[[i]]
  }
  cat("\n")
  
  aeics=aeics[!sapply(aeics,is.null)]
  
  if(conv2df)
    aeics=lapply(aeics,function(y) do.call("rbind",lapply(names(y),function(x) cbind(Sid=x,data.frame(y[[x]])))))
  
  return(aeics)
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Merge ROIs/EICs from several files not very fast plus pbs with duplicated roi at the end  and/or scan
#' 
#' @param files list of files
#' @param listRois ROI/EICs names to recoup across files
#' @param parDeco par deco
#' @param lsamp2use sample to use to compute stats
#' @param minHeight minimum peak height in lsamp2use
#' @param addrt extra padding around frames
#' @param doPlot should it be plotted
#' 
#' @import matrixStats
#' @import MASS
#' @export

eicGather2<-function(fileres,listRois,parDeco,lsamp2use,minHeight=parDeco$minHeightMS1,addrt=NULL,doPlot=T){
  
  
  .infctreframe<-function(eic,lsamp2use,parDeco,minHeight,addrt,doPlot,main=""){
#    for(ix in names(leics)){
 #     eic=leics[[ix]]
      l=which(eic$Sid%in%lsamp2use & !is.na(eic$mz))
      ieic=eic[l,]
      rmz=range(ieic[,"mz"])
      mindmz=min(c(mean(rmz)*parDeco$ppm*10^-6,parDeco$dmz))/2
      
      if(doPlot){
        ymz=range(pretty(range(eic$mz,na.rm=TRUE)+c(-0.001,0.001)))
        yrt=range(pretty(range(eic$nrt,na.rm=TRUE)+c(-0.2,0.2)))
        l=which(!is.na(eic$mz))
        k <- MASS::kde2d(eic$mz[l], eic$nrt[l],h = c(0.001,0.05),lims = c(ymz,yrt))
        #      k$z=k$z #*sum(!is.na(eic$mz))/sum(k$z)
        k$z[k$z<0.001]=NA
        image(k, col=colorRampPalette(rev(brewer.pal(9,'Greys')))(11)[11:6],main=main)
        abline(v=seq(min(ymz),max(ymz),median(ymz)*5*10^-6))
        ic=as.numeric(cut(ieic$y,c(0,parDeco$minHeightMS1,minHeight,Inf)))
        l=which(ieic$y>=parDeco$minHeightMS1)
        points(ieic$mz[l],ieic$nrt[l],pch=16,col=ic[l],cex=c(.2,1,.8)[ic[l]])
        #      abline(v=range(ieic$mz,na.rm=TRUE))
      }
      lusid=unique(ieic$Sid[ieic$y>=parDeco$minHeightMS1])
      ### loop over each sample
      ncons=ifelse(is.null(parDeco$ncons),3,parDeco$ncons)
      ppmcons=ifelse(is.null(parDeco$ppm2),parDeco$ppm,parDeco$ppm2)
      minNoise=ifelse(is.null(parDeco$minHeightMS1cons),parDeco$minHeightMS1,parDeco$minHeightMS1cons)
      allre=NULL
      for(ii in lusid){
        l2=which(ieic$Sid== ii& ieic$y>parDeco$minHeightMS1)
        llx<-.MGinquickSplit(ieic[l2,],ncons=ncons,dppm =ppmcons,dmz=parDeco$dmz,minNoise=minNoise)[[1]]
        if(length(llx)==0) next
        ldrmz=sapply(llx,function(x) diff(range(ieic[l2[x],"mz"])))
        l=which(ldrmz<=max(quantile(ldrmz,.75),mindmz))
        mllx=lapply(GRMeta:::.GRmergellx(llx[l]),function(x) l2[unique(x)])
        re=do.call("rbind",lapply(mllx,function(x) c(range(ieic$mz[x]),range(ieic$nrt[x]),max(ieic$y[x]),ieic$mz[x[which.max(ieic$y[x])]])))
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
      newwin=newwin[which(newwin[,'int']>=minHeight),,drop=F]
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
        l=which(eic$nrt>=(x[3]-addrt) & eic$nrt<=(x[4]+addrt))
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
      
      if(doPlot) apply(newwin,1,function(x) rect(x[6],x[3]-1*addrt,x[7],x[4]+1*addrt,border = 4,lwd=2))
  
    return(leics)
  }
  
  aeics=aeics0=eicGather(fileres,listRois,conv2df = T)
  
  # (ix=names(aeics)[17])
  if(is.null(addrt))  addrt=2*parDeco$psdrt*parDeco$span+parDeco$drt*2 ## add 0.5 minutes for bsl corr stuff
  aeics=aeics0
  ## 1) first pass screening set of eics first
  llneic=list()
  for(ix in names(aeics)){
    eic=aeics[[ix]]
    l=which(eic$Sid%in%lsamp2use)
    if(max(eic$y[l],na.rm=T)<minHeight){
      aeics[[ix]]=NULL
      next
    }
    l=which(eic$Sid%in%lsamp2use & eic$y>=parDeco$minHeightMS1)
    ll=.GRsplist(eic$nrt[l],l,d = 2*addrt)
    ll=ll[which(sapply(ll,function(x) max(eic$y[x]))>=minHeight)]
    ll=lapply(ll,function(x) range(eic$nrt[x])+c(-1,1)*addrt)
    if(length(ll)>1) print(c("",ix))
    neic=lapply(ll,function(x) eic[which(eic$nrt>=x[1] & eic$nrt<=x[2]),,drop=F])
    if(length(neic)>1) names(neic)=paste0(ix,"-",1:length(neic)) else names(neic)=ix
    aeics=c(aeics,neic)
    aeics[[ix]]=NULL
  }
  if(length(aeics)==0) return(NULL)
  
  ## 2) second pass screening set of eics first
  l2chk=c()
  neics=list()
  for(ix in names(aeics)){
    cat(".")
    eic=aeics[[ix]]
    l=which(eic$Sid%in%lsamp2use)
    
    rmz=range(eic$mz[l],na.rm=T)
    rmz0=range(eic$mz,na.rm=T)
    rmz2=matrixStats:::weightedMedian(eic$mz[l],eic$y[l],na.rm=T)
    rmz2=range(rmz2*(1+c(-1,1)*parDeco$ppm*10^-6),rmz2+c(-1,1)*parDeco$dmz)
    isok=max(table(eic$scan[l],eic$Sid[l])) ## no duplicates
    isok1=(rmz[1]>rmz2[1] & rmz[2]<rmz2[2]) ## all sample are within parDeco
    isok2=(rmz0[1]>rmz2[1] & rmz0[2]<rmz2[2]) ## all datapoints are within the sample range
    if(isok==1 & isok2 & isok1){ ## the 3 
      l=which(eic$Sid%in%lsamp2use & !is.na(eic$mz))
      ieic=eic[l,]
      lusid=names(sort(table(ieic$Sid[ieic$y>=parDeco$minHeightMS1]),dec=T))
      cond=TRUE
      while(cond & length(lusid)>0){
        ii=lusid[1]
        l2=which(ieic$Sid== ii& ieic$y>=parDeco$minHeightMS1)
        llx<-.MGinquickSplit(ieic[l2,],ncons=3,dppm = parDeco$ppm,dmz=parDeco$dmz,minNoise=parDeco$minHeightMS1)[[1]]
        if(length(llx)>0) cond=FALSE ## at least one sample contains 3 consecutive scans
        lusid=lusid[-1]
      }
      if(cond){
        print(ix)
        aeics[[ix]]=NULL
      }
      next
    }
    re=.infctreframe(eic,lsamp2use,parDeco,minHeight,addrt,doPlot,main = ix)
    if(!is.null(re)){
      names(re)=paste0(ix,"-",1:length(re))
      neics=c(neics,re)
    }
    aeics[[ix]]=NULL
    
  }
  if(length(neics)>0) aeics=c(aeics,neics)
  if(length(aeics)==0) return(NULL)
  
  # ########### 3) reduce the m/z and rt range if too wide
  # if(length(l2chk)>0){
  #   leics=.infctreframe(aeics[l2chk],lsamp2use,parDeco,minHeight,addrt,doPlot)
  #   aeics=aeics[!names(aeics)%in%l2chk]
  # if(length(leics)>0) aeics=c(aeics,leics)
  # }
  
  ########## check if rois overlap
  eicmat2=data.frame(do.call("rbind",lapply(aeics,function(x) c(round(range(x$nrt),4),round(range(x$mz,na.rm=T),6),
                                                                round(matrixStats:::weightedMedian(x$mz,x$y,na.rm=T),6)))))
  names(eicmat2)=c("rtmin","rtmax","mzmin","mzmax","mz50")
  eicmat2$RoiId=sprintf("R%.5f@%.2f-%.2f",eicmat2$mz50,eicmat2$rtmin,eicmat2$rtmax)
  llover=.GRisover2(eicmat2$mzmin,eicmat2$mzmax,eicmat2$rtmin,eicmat2$rtmax,T,10^-4,parDeco$psdrt)
  
  l2reframe=which(sapply(llover,length)>1)
  if(length(l2reframe)>0){
    for(ix in l2reframe){
      leics=do.call("rbind",aeics[llover[[ix]]])
      leics=leics[which(!duplicated(leics)),]
      leics=leics[order(leics$Sid,leics$nrt),]
      re=.infctreframe(leics,lsamp2use,parDeco,minHeight,addrt,doPlot)
      if(!is.null(re)) aeics=c(aeics,re)
    }
    aeics=aeics[-unlist(llover[l2reframe])]
  }
  
  eicmat3=data.frame(do.call("rbind",lapply(aeics,function(x) c(round(range(x$nrt),4),round(range(x$mz,na.rm=T),6),
                                                                round(matrixStats:::weightedMedian(x$mz,x$y,na.rm=T),6)))))
  names(eicmat3)=c("rtmin","rtmax","mzmin","mzmax","mz50")
  eicmat3$RoiId=sprintf("R%.5f@%.2f-%.2f",eicmat3$mz50,eicmat3$rtmin,eicmat3$rtmax)
  #llover=.GRisover2(eicmat3$mzmin,eicmat3$mzmax,eicmat3$rtmin,eicmat3$rtmax,T,10^-4,parDeco$psdrt)
  ldups=which(!duplicated(eicmat3))
  eicmat3=eicmat3[ldups,]
  aeics=aeics[ldups]
  eicmat3$RoiIdo=names(aeics)
  rownames(eicmat3)=names(aeics)=eicmat3$RoiId
    
  # if(length(llneic)>0) aeics=c(aeics,llneic)
  if(length(aeics)==0) return(NULL)
  
  lso=order(eicmat3$mz50,eicmat3$rtmin)
  aeics=aeics[lso]
  eicmat3=eicmat3[lso,]
  list(Eic=aeics,EicDef=eicmat3)
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' fill holes etc...
#' 
#' @param xeic EIC
#' @param drt pseudo delta rt
#' @param span span
#' @param bw bandwith calculation
#' @param minNoise noise intensity
#' @param retbsl should baseline be computed
#' @keywords internal
#' 
#' @export
.MGfill<-function(xeic, drt, span = 5, bw = drt * span * 2, minNoise = 5000,retbsl=TRUE,minpts=3){
  
  # span = 5;bw = drt * span * 2;minNoise = 10000;minHeight = NA
  
  lnnas = which(!is.na(xeic[, "y"]) & xeic[, "y"] >= minNoise)
  if(length(lnnas)<minpts) return(NULL)
  segpks = GRMeta:::.GRsplist(xeic[lnnas, "rt"], lnnas, d = span *drt/1.99, ismass = F)
  xrt = t(sapply(segpks, function(y) range(xeic[y, "rt"])))
  xrt[, 1] = xrt[, 1] - bw/2
  xrt[, 2] = xrt[, 2] + bw/2
  if (nrow(xrt) > 1) {
    ll2merge = GRMeta:::.GRisover(xrt[, 1], xrt[, 2], T, bw - drt)
    xrt = do.call("rbind", lapply(ll2merge, function(x) range(xrt[x,])))
  }
  xeic = cbind(xeic, isin = 0)
  for (i in 1:nrow(xrt)) xeic[xeic[, "rt"] >= xrt[i, 1] & xeic[, "rt"] <= xrt[i, 2], "isin"] = 1
  xeic = cbind(xeic, y2 = xeic[, "y"])
  xeic[which(xeic[, "isin"] == 0 | xeic[, "y"] <= minNoise), "y2"] = minNoise
  if (is.na(xeic[1, "y2"])) xeic[1, "y2"] = minNoise
  if (is.na(xeic[nrow(xeic), "y2"]))  xeic[nrow(xeic), "y2"] = minNoise
  xeic = cbind(xeic, toimp = 0)
  lnna = which(is.na(xeic[, "y2"]))
  if (length(lnna)) {
    l2rep = GRMeta:::.GRsplist(lnna, lnna, d = 1.1)
    for (y in l2rep) {
      df = data.frame(Y = log(xeic[range(y) + c(-1, 1), "y2"]), RT = xeic[range(y) + c(-1, 1), "rt"])
      xeic[y, "y2"] = exp(predict(lm(Y ~ RT, df), newdata = data.frame(RT = xeic[y, "rt"])))
      xeic[y, "toimp"] = 1
    }
  }
  
  newx = .MGdoksmooth(x = xeic[, "rt"], y = xeic[, "y2"], missc = xeic[, "toimp"] == 1, bw = bw, drt = NA)
  lsc=min(xeic[,"nscan"]):max(xeic[,"nscan"])
  #newx=data.frame(nscan=lsc,rt=lsc*drt,y2=approx(newx$x,newx$y,lsc*drt)$y)  ### NAs at the end-> pb of overlap iwth new vs old scan rt!!!
  newx=data.frame(xeic[match(lsc,xeic[,"nscan"]),c("scan","mz","y","id")],nscan=lsc,rt=lsc*drt,y2=approx(newx$x,newx$y,lsc*drt)$y)
  newx=newx[!is.na(newx$y2),]
  
  if(!retbsl) return(newx)
  y=newx$y2
  n2pad=131
  y=c(rep(minNoise,n2pad-span*3-1),rep(y[1],span*3+1),y,rep(rev(y)[1],span*3+1),rep(minNoise,n2pad-span*3-1))
  bsl=GRMeta:::.GRbslrf(1:length(y),y,NoXP = NULL)
  bsl$fit[bsl$fit<minNoise]=minNoise
  bslscore <- (y - bsl$fit)/max(bsl$sigma, 10^-3)
  bslscore[which(abs(bslscore) > 10)] = sign(bslscore[which(abs(bslscore) > 10)]) * 10
  
  newx$bsl = bsl$fit[n2pad+(1:length(newx$y2))]
  newx$bslc = bslscore[n2pad+(1:length(newx$y2))]
  newx
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Fetch ROIs/EICs across files
#'
#' @param listEics list of eics to be processed
#' @param parDeco deconvolution paramters
#' @param refFct cleaning function
#' @param nSlaves number of slaves for parallel
#' @param outfile name of the outfile
#' @param what columns to export
#' @param chunk chunk size for parallel processing
#' @param verbose verbose
#' @import foreach
#' @import doParallel
#' @return list of eics
#' @export
.MGfillMulti<-function(listEics,parDeco,refFct=metaboGoS:::.MGfill,nSlaves=1,outfile=NULL,what=c("y","y2","mz"),chunk=50,verbose=TRUE){
  
  allxeic=list()
  ll=names(listEics)#eicmat$RoiId#[101:150]
  ll=split(ll, ceiling(seq_along(ll)/chunk))
  n0=length(listEics)
  listEics=lapply(ll,function(x) listEics[x])
  cat("Processing ", n0, " Roi/EICs split in ",length(ll)," chunk of avg. size=",round(mean(sapply(ll,length)),1),sep="")
  # 
  # 
  # 
  # 
  #   lperc=""
  # if(length(ll)>20) lperc=ll[round(seq(1,length(ll),length=12)[2:11])]
  # if(length(ll)>50) lperc=ll[round(seq(1,length(ll),length=22)[2:21])]
  # 
  if(nSlaves>1)   nSlaves=max(1, min(nSlaves,detectCores()-1))
  if(nSlaves>1){
    if(nSlaves>length(listEics)) nSlaves=length(listEics)
    clProc<-makeCluster(nSlaves)
    doParallel:::registerDoParallel(clProc)
    cat(" -- registering ",nSlaves," clusters\n",sep="")
  }
  
  .infctMGfillMulti<-function(leics,parDeco,what,tempout=TRUE,verbose=TRUE){
    
    eiclist=list()
    for(iroi in names(leics)){ 
      if(verbose) cat(".")
      ieic=list()
      xeic=leics[[iroi]]
      for(ii in unique(xeic$Sid))
        ieic[[ii]]=.MGfill(xeic[xeic$Sid==ii,],drt=parDeco$psdrt,span = parDeco$span,minNoise=parDeco$minNoiseMS1,bw=parDeco$bw,retbsl = FALSE,minpts = 3)[,c(what,"nscan",'rt')]
      if(length(ieic)==0) next
      lsc=max(sapply(ieic,function(x) min(x[,"nscan"]))):min(sapply(ieic,function(x) max(x[,"nscan"])))
      ieic=lapply(ieic,function(x) x[match(lsc,x[,"nscan"]),what])
      ieic=lapply(what,function(iwhat) t(sapply(ieic,function(x) x[,iwhat]))) 
      ieic=lapply(ieic,function(x){colnames(x)=lsc;x})
      names(ieic)=what
      eiclist[[iroi]]=ieic
      rm(list=c('xeic','ieic'))
    }
    if(tempout){
      if(length(eiclist)==0) return(NULL)
      outfile=paste0(tempdir(),"/tmp",gsub("[A-Z@-\\.]","",names(eiclist)[1]),".rda")
      save(file=outfile,eiclist)
      return(outfile)
    }
    invisible(eiclist)
  }
  
  
  ### Parallele bit
  if(nSlaves>1){
    #llre=foreach(i = ll, .export = fct2exp,.packages = c("igraph","xcms","GRMeta"), .verbose =F)  %dopar%{
    ltmpfiles=foreach(ieics = listEics,.export = ".infctMGfillMulti",.packages = c("metaboGoS"), .verbose =F)  %dopar%{
      .infctMGfillMulti(ieics,parDeco,what,tempout=TRUE,verbose = FALSE)
    }
    ltmpfiles=unlist(ltmpfiles)
  }
  
  ## Serial bit
  if(nSlaves<=1){
    ltmpfiles=c()
    for(ii in 1:length(listEics)){
      if(verbose) cat("Chunk ",ii,": ",sep="")
      ltmpfiles[ii]=.infctMGfillMulti(listEics[[ii]],parDeco,what,tempout=TRUE,verbose = verbose)
    }
  }
  
    ltmpfiles=ltmpfiles[file.exists(ltmpfiles)]
  ### sort the rest out
  allxeic=list()
  for(ifile in ltmpfiles){
    load(ifile)
    if(!exists("eiclist")) next
    if(!is.list(eiclist)) next
    allxeic=c(allxeic,eiclist)
    rm(list="eiclist")
  }
  
  cat('---> Found ',length(allxeic),"/",length(listEics),"\n",sep="")
  if(nSlaves>1) stopCluster(clProc)
  # if(save & is.null(outfile)) save(file=paste0(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", infile),".rda"),allxeic,sc2nrt)
  if(!is.null(outfile)) save(file=outfile,allxeic)
  invisible(allxeic)
  
}


