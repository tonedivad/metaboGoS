### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Ridiculusly complicated mergingof ROIS
#' 
#' @param matroi Data frame of ROIs
#' @param currstats Stats 
#' @param dppm dppm
#' @param dmz dmz
#' @param lbw band withs
#' @param maxIter maxIter
#' 
#' 
#' @export
mergeROIs<-function(matroi,dppm=2.5,dmz=0.001,lbw=seq(dppm*2,dppm/2,length.out = 5),maxIter=9,
                    indicatorVec=c(sid="Sid",grp="GrpEIC",height='intensity',
                                   frmin="mz10",frcen="mz50",frmax="mz90",
                                   mzmin="mzmin",mz="mz50",mzmax="mzmax",
                                   rtmin="rtmin",rtmax="rtmax",rt="rt")){
  
  # indicatorVec=c(sid="Sid",grp="GrpEIC",height='intensity',frmin="mz10",frcen="mz50",frmax="mz90",mzmin="mzmin",mz="mz50",mzmax="mzmax",rtmin="rtmin",rtmax="rtmax",rt="rt")
  # dppm=2.5;dmz=0.001;lbw=seq(dppm*2,dppm/2,length.out = 5);maxIter=9
  
  ## rough grouping on mz
  currsplit <- .GRsplistMZ(matroi[,indicatorVec["frcen"]],  dppm = NA,dmz=20*dmz)
  if(any(sapply(currsplit,is.list))) currsplit=unlist(currsplit,recursive = FALSE)
  tmp=do.call("rbind",lapply(1:length(currsplit),function(x) cbind(x,currsplit[[x]])))
  matroi[tmp[,2],indicatorVec["grp"]]=tmp[,1]
  currois0=matroi
  currstats0=infctMRoi.compstats(matroi,indicatorVec)#,mzmin = "mzmin",mzmax = "mzmax")
  
  ## make sure there is no overlap for the frame min/max
  re=infctMRoi.mergemzrt(currois0,currstats0,typ = indicatorVec[c("frmin","frmax")],indicatorVec,dppm = dppm,dmz=dmz/2)
  if(!is.null(re)){currois0=re$rois;currstats0=re$st}
  cat(" --> step1: ",nrow(currstats0)," rois out of ",nrow(currois0),"\n",sep="")
  
  ## step 2 merge close mz/overlap
  crois=currois0
  cstats=currstats0
  doLoop=T;iter=0
  cat(" --> step2: ",sep="")
  while(doLoop & iter<maxIter){
    iter=iter+1
    cat(iter,sep="")
    re1=infctMRoi.mergeclosemz(crois,cstats,indicatorVec,thrMZ=dmz/2,dmz=dmz,dppm=dppm)
    if(!is.null(re1)){
      crois=re1$rois[!duplicated(re1$rois),];cstats=re1$stats
      cat("-",nrow(crois)," ",sep="")
    } else doLoop=F
  }
  cat("\n --> ",nrow(crois)," rois / ",nrow(cstats)," eic groups\n",sep="")
  
  ## step 3 further split based on mz and RT
  re2=infctMRoi.splitmzrt(crois,cstats,indicatorVec)
  if(!is.null(re2)){crois=re2$rois;cstats=re2$stats}
  re1=infctMRoi.mergeclosemz(crois,cstats,indicatorVec,thrMZ=dmz,dmz=dmz,dppm=dppm) ## just in case
  if(!is.null(re1)){crois=re1$rois;cstats=re1$stats}
  ## step 3b make sure mz are ordered and grpEIC from 1 to max
  cstats=cstats[order(cstats[,indicatorVec["frcen"]]),]
  crois[,indicatorVec["grp"]]=as.numeric(factor(crois[,indicatorVec["grp"]],unique(cstats[,indicatorVec["grp"]])))
  cstats=infctMRoi.compstats(crois,indicatorVec)
  cat(" --> step3: ",nrow(crois)," rois / ",nrow(cstats)," eic groups\n",sep="")
  
  ## step 4 refine based on density
  cat(" --> step4: ",sep="")
  for(ibw in lbw){
    cat(" ",ibw,"-",sep="")
    l2split=cstats$GrpEIC[which((cstats$dppm>=(2*dppm) | cstats$dppm0>=(2*dppm)) | cstats$mndups>1)]
    if(length(l2split)==0) next
    re3=infctMRoi.splitmzdens(crois,cstats,lurois = l2split,indicatorVec,bw=ibw)
    if(!is.null(re3)){crois=re3$rois;cstats=re3$stats;cat(nrow(cstats),sep="")}
  }
  cstats=cstats[order(cstats[,indicatorVec["frcen"]]),]
  crois[,indicatorVec["grp"]]=as.numeric(factor(crois[,indicatorVec["grp"]],unique(cstats[,indicatorVec["grp"]])))
  cstats=infctMRoi.compstats(crois,indicatorVec)
  cat(" --> step4: ",nrow(crois)," rois / ",nrow(cstats)," eic groups\n",sep="")
  
  eicmat=data.frame(RoiId=NA,ovl=NA,
                    rt=round(floor(cstats[,indicatorVec['rt']]*10000)/10000,4),
                    rtmin=round(floor(cstats[,indicatorVec['rtmin']]*10000)/10000,4),
                    rtmax=round(ceiling(cstats[,indicatorVec['rtmax']]*10000)/10000,4),
                    mz=round((cstats[,indicatorVec['mz']]+cstats[,indicatorVec['frcen']])/2,6),
                    mz10=round(cstats[,indicatorVec['frmin']],6),
                    mz90=round(cstats[,indicatorVec['frmax']],6),mzmin=NA,mzmax=NA)
  mmz=do.call("rbind",lapply(1:nrow(eicmat),function(i) range(eicmat[i,c("mz","mz","mz10","mz90")]+c(-1,1,-.5,.5)*dppm*10^-6*eicmat$mz[i])))
  eicmat$mzmin=round(mmz[,1],6)
  eicmat$mzmax=round(mmz[,2],6)
  eicmat=eicmat[order(eicmat$mz),]
  llsplit=.GRisover2(eicmat$mzmin,eicmat$mzmax,eicmat$rtmin,eicmat$rtmax,ret=TRUE)
  llsplit=do.call("rbind",lapply(1:length(llsplit),function(x) cbind(llsplit[[x]],x)))
  eicmat$ovl[llsplit[,1]]=llsplit[,2]
  eicmat$RoiId=sprintf("R%.5f@%.2f-%.2f",eicmat$mz,eicmat$rtmin,eicmat$rtmax)
  invisible(eicmat)
 # invisible(list(eicmat=eicmat,crois=crois))
}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Gather the rois for each resfile
#' 
#' @param fileres List of files
#' @param what what
#' 
#' @keywords internal
#' 
#' @export
infctMRoi.get<-function(fileres,what='RoiInfos'){
  
  fileres=fileres[file.exists(fileres)]
  if(length(fileres)==0) stop('No file found')
  if(is.null(names(fileres))) names(fileres)=paste0("S",1:length(fileres))
  
arois=list()
for(isid in names(fileres)){
  cat(".")
  load(fileres[isid])
  y=obj[[what]]
  if(is.null(y)) next
  y$Sid=isid
  arois[[isid]]=y
}
arois=do.call("rbind",arois)
if(any(grepl("mz",names(arois)))) arois=arois[order(rowMeans(arois[,grep("mz",names(arois))])),]
rownames(arois)=NULL
return(arois)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Wrapper for infctMRoi.chkclosemz: loop over GrpEIC
#' 
#' @param currrois Data frame of ROIs
#' @param currstats Stats 
#' @param thrMZ overlap in mz
#' @param dmz dmz
#' @param dppm dppm
#' 
#' @keywords internal
#' 
#' @export
infctMRoi.mergeclosemz<-function(currrois,currstats,indicatorVec,thrMZ=0.001,dmz=0.0005,dppm=2.5){
  
  l=which(currrois[,indicatorVec['grp']]%in%currstats[,indicatorVec['grp']])
  ll2chk=tapply(l,currrois[l,indicatorVec['grp']],c)
  l2merge=list()
  for(lx in ll2chk){
    tmproi=currrois[lx,]
    l2m=infctMRoi.chkclosemz(tmproi,indicatorVec,thrMZ=thrMZ,dmz=dmz,dppm=dppm)
    if(length(l2m)==0) next
    l2merge=c(l2merge,lapply(l2m,function(y) lx[y]))
  }
  #print(l2merge)
  if(length(l2merge)==0) return(NULL)
  narois=do.call("rbind",lapply(l2merge,function(x) infctMRoi.merge(currrois[x,])))
  narois=rbind(currrois[-unlist(l2merge),],narois)
  narois=narois[order(narois[,indicatorVec['grp']],narois[,indicatorVec['sid']],narois[,indicatorVec['height']]),]
  narois[,indicatorVec['grp']]=as.numeric(factor(narois$GrpEIC))
  nstats=infctMRoi.compstats(narois,indicatorVec)
  invisible(list(rois=narois,stats=nstats))
}

infctMRoi.merge<-function(idf,orderby='intensity',
                      lminvar=c("rtmin","mz10","mzmin"),
                      lmaxvar=c("rtmax","mz90","mzmax","coda"),
                      lmedvar=c("mz50")){
  idf=idf[order(idf[,orderby],decreasing = TRUE),]
  for(i in lminvar[lminvar%in%names(idf)]) idf[1,i]=min(idf[,i])
  for(i in lmaxvar[lmaxvar%in%names(idf)]) idf[1,i]=max(idf[,i])
  for(i in lmedvar[lmedvar%in%names(idf)]) idf[1,i]=median(idf[,i])
  idf[1,]
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' 1 Merge based on dmz or overlap in mz 2 merge on mz/rt  wrapper for infctMRoi.loopmzrt
#' 
#' @param currrois Data frame of ROIs
#' @param currstats Sample id indicator 
#' @param typ mz column indicator
#' @param mzmin mzmin
#' @param mzmax mzmax
#' @param thrMZ overlap in mz
#' @param dmz dmz
#' @param dppm dppm
#' @return list of entries to be merged in tmproi - may be null
#' @keywords internal
#' 
#' @export

infctMRoi.mergemzrt<-function(currrois,currstats,typ="mzmed",indicatorVec,dppm=2.5,dmz=0.0005,thrMZ=0.001,thrRT=0,maxiter=6){
  
  if(length(typ)==1) csplit<-.GRsplistMZ(currstats[,typ[1]],dppm=dppm,dmz = dmz)
  if(length(typ)==2) csplit<-.GRisover(currstats[,typ[1]],currstats[,typ[2]],ret=TRUE,thr = thrMZ)
  
  grp=indicatorVec["grp"]
  n=sapply(csplit,length)
  l2split=which(n>1)
  nsplit=list()
  cid=max(currstats[,grp])
  for(idx in l2split){
    lx=which(currrois[,grp]%in%currstats[csplit[[idx]],grp])
    tmproi=currrois[lx,]
    v=infctMRoi.loopmzrt(tmproi,indicatorVec,thrMZ=thrMZ,thrRT=thrRT,maxiter=maxiter)
    newv=unlist(lapply(1:length(v),function(x) rep(x,length(v[[x]]))))
    if(sum(diag(table(newv,tmproi[,grp])))==length(lx)) next
    nsplit=c(nsplit,list(cbind(x=lx,oGrpEIC=tmproi[,grp],nGrpEIC=newv+cid)))
    cid=cid+max(newv)
  }
  if(length(nsplit)==0) return(NULL)
  nsplit=do.call("rbind",nsplit)
  currrois[nsplit[,1],grp]=nsplit[,3]
  narois=currrois[order(round(currrois[,indicatorVec['frcen']]),currrois[,grp],currrois[,indicatorVec['sid']],currrois[,indicatorVec['height']]),]
  narois[,grp]=as.numeric(factor(narois[,grp]))
  nstats=infctMRoi.compstats(narois,indicatorVec)
  invisible(list(rois=narois,stats=nstats))
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Double merge on mz/rt - wrapper for infctMRoi.loopmzrt
#' 
#' @param currrois Data frame of ROIs
#' @param currstats Sample id indicator 
#' @param indicatorVec indicator
#' @param lurois vector of rois
#' @param cid grp max value to start from
#' @param thrMZ overlap in mz
#' @param thrRT overlap in rt
#' @return list of entries to be merged in tmproi - may be null
#' @keywords internal
#' 
#' @export

infctMRoi.splitmzrt<-function(crois,cstats,indicatorVec,lurois=unique(cstats[,indicatorVec['grp']]),cid=max(cstats[,indicatorVec["grp"]]),thrMZ=0.001,thrRT=0.001){
  nsplit=list()
  grp=indicatorVec["grp"]
  for(ieic in lurois){
    lx=which(crois[,grp]==ieic)
    tmproi=crois[lx,]
    v=.GRisover2(tmproi[,indicatorVec["frmin"]],tmproi[,indicatorVec["frmax"]],tmproi[,indicatorVec["rtmin"]],tmproi[,indicatorVec["rtmax"]],
                 retOne = TRUE,thr1 = thrMZ,thr2=thrRT)
    if(length(v)<2) next
    newv=do.call("rbind",lapply(1:length(v),function(x) cbind(x=lx[v[[x]]],oGrpEIC=tmproi[v[[x]],grp],nGrpEIC=x+cid)))
    nsplit=c(nsplit,list(newv))
    cid=max(newv[,3])
  }
  if(length(nsplit)==0) return(NULL)
  nsplit=do.call("rbind",nsplit)
 # l2rm=nsplit[,1]
  nrois=crois[nsplit[,1],]
  nrois[,grp]=nsplit[,3]
  #nrois$oGrpEIC=nsplit[,2]
  nstats=infctMRoi.compstats(nrois,indicatorVec)
  luold=unique(nsplit[,2])
  ncstats=rbind(cstats[!cstats[,grp]%in%luold,],nstats)
  ncrois=rbind(crois[!crois[,grp]%in%luold,],nrois)
  return(list(rois=ncrois,stats=ncstats))
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Split ROI based on density mz  
#' 
#' @param crois ROIs
#' @param cstats Stats  
#' @param lurois vector of grp ids be splitted
#' @param indicatorVec indicatorVec
#' @param bw ppm bandwidth parameter
#' @param cid Group EIC numeric id
#' @return list ROI/stats after splitting crois -> replace crois/cstats
#' @keywords internal
#' 
#' @export

infctMRoi.splitmzdens<-function(crois,cstats,lurois,indicatorVec,bw=5,cid=max(cstats[,indicatorVec['grp']])){
  nsplit=list()
  grp=indicatorVec['grp']
  mz=indicatorVec['mz']
  for(ieic in lurois){
    lx=which(crois[,grp]==ieic)
    tmproi=crois[lx,]
    v=.GRsplitOne(tmproi[,mz],bw = bw,ismass = TRUE)$sp
    if(length(v)<2) next
    newv=do.call("rbind",lapply(1:length(v),function(x) cbind(x=lx[v[[x]]],oGrpEIC=tmproi[v[[x]],grp],nGrpEIC=x+cid)))
    nsplit=c(nsplit,list(newv))
    cid=max(newv[,3])
  }
  if(length(nsplit)==0) return(NULL)
  nsplit=do.call("rbind",nsplit)
  # l2rm=nsplit[,1]
  nrois=crois[nsplit[,1],]
  nrois[,grp]=nsplit[,3]
  #nrois$oGrpEIC=nsplit[,2]
  nstats=infctMRoi.compstats(nrois,indicatorVec)
  luold=unique(nsplit[,2])
  ncstats=rbind(cstats[!cstats[,grp]%in%luold,],nstats)
  ncrois=rbind(crois[!crois[,grp]%in%luold,],nrois)
  return(list(rois=ncrois,stats=ncstats))
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Merge slices from the same sample within a single ROI 
#' 
#' @param tmproi Data frame of ROIs
#' @param sid Sample id indicator 
#' @param mzmed mzmed
#' @param mzmin mzmin
#' @param mzmax mzmax
#' @param thrMZ overlap in mz
#' @param dmz dmz
#' @param dppm dppm
#' @return list of entries to be merged in tmproi - may be null
#' @keywords internal
#' 
#' @export
infctMRoi.chkclosemz<-function(tmproi,indicatorVec,thrMZ=0.001,dmz=0.0005,dppm=2.5){
  
  sid=indicatorVec["sid"]
  mzmed=indicatorVec["frcen"]
  mzmin=indicatorVec["frmin"]
  mzmax=indicatorVec["frmax"]
  l2merge=list()
  tab=table(tmproi[,sid])
  # print(tab[tab>1])
  lsid=names(which(tab>1))
  for(isid in lsid){
    l=which(tmproi[,sid]==isid)
    itmproi=tmproi[tmproi[,sid]==isid,]
    ill=lapply(.GRsplistMZ(itmproi[l,mzmed],dppm = dppm,dmz = dmz),function(x) l[x])
    ill=ill[lapply(ill,length)>1]
    if(length(ill)==0) next
    ill=unlist(lapply(ill,function(x) .GRisover(tmproi[x,mzmin],tmproi[x,mzmax],retOne = TRUE,thr = thrMZ)),recursive = FALSE)
    ill=lapply(ill,function(x) l[x])
    l2merge=c(l2merge,ill[lapply(ill,length)>1])
  }
  return(l2merge)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Iterative merging of slices with overlaping mz/rt within a single ROI regardless sample Id
#' 
#' @param tmproi Data frame of ROIs
#' @param indicatorVec column indicator vector
#' @param thrMZ overlap in mz
#' @param thrRT overlap in rt
#' @param maxiter maximum iterations, usually done within 1-2
#' @return list of entries in tmproi
#' @keywords internal
#' 
#' @export
infctMRoi.loopmzrt<-function(tmproi,indicatorVec,thrMZ=0.001,thrRT=0.0001,maxiter=6){
  v=list(1:nrow(tmproi))
  rtmin=indicatorVec['rtmin']
  rtmax=indicatorVec['rtmax']
  mzmin=indicatorVec['frmin']
  mzmax=indicatorVec['frmax']
  
  nv=1;doLoop=TRUE;iter=0
  while(doLoop & iter<=maxiter){
    iter=iter+1
    v=.GRisover(tmproi[,mzmin],tmproi[,mzmax],retOne = TRUE,thr = thrMZ)
    v=unlist(lapply(v,function(y){
      z=.GRisover(tmproi[y,rtmin],tmproi[y,rtmax],retOne = TRUE,thr = thrRT)
      lapply(z,function(i) y[i])
    }),rec=F,use=F)
    doLoop=length(v)!=nv
    nv=length(v)
  }
  return(v)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' Compute stats for each group of EIC
#' 
#' @param arois Data frame with ROI definition for each sample
#' @param grp EIC group column indicator 
#' @param sid Sample id column indicator 
#' @param mz average mz column indicator
#' @param mzmin mzmin column indicator 
#' @param mzmax mzmax column indicator 
#' @param rtmin rtmin column indicator 
#' @param rtmax rtmax column indicator 
#' @return data.frame with summary of each group of EIC
#' @keywords internal
#' 
#' @export
infctMRoi.compstats<-function(arois,indicatorVec){
  # grp="GrpEIC";mzmed="mz50";mzmin="mz10";mzmax="mz90";rtmin="rtmin";rtmax="rtmax"
  # currsplit=tapply(1:nrow(arois),arois[,grp],c)
  # currstats=t(sapply(currsplit,function(x){
  #   v=10^6*(range(arois$mz50[x])/median(arois$mz50[x])-1)
  #   vrt=range(arois$rtmin[x],arois$rtmax[x])
  #   c(median(arois$mz50[x]),diff(v),v,median(arois$rt[x]),vrt,length(unique(arois$Sid[x])),max(table(arois$Sid[x])),length(x))
  # }))
  sid=indicatorVec['sid']
  grp=indicatorVec['grp']
  mz=indicatorVec['mz']
  rt=indicatorVec['rt']
  rtmin=indicatorVec['rtmin']
  rtmax=indicatorVec['rtmax']
  mzmed=indicatorVec['frcen']
  mzmin=indicatorVec['frmin']
  mzmax=indicatorVec['frmax']
  
  currstats=data.frame(do.call("rbind",tapply(1:nrow(arois),arois[,grp],function(x){
 #   v0=10^6*(range(arois[x,mzmed])/median(arois[x,mz])-1)
    v=10^6*(range(arois[x,mzmed])/median(arois[x,mzmed])-1)
    vrt=range(arois[x,rtmin],arois[x,rtmax])
    c(arois[x[1],grp],
   #  median(arois[x,mzmed]),diff(v0),
      median(arois[x,mzmed]),diff(v),v,
      range(arois[x,mzmin],arois[x,mzmax]),
      median(arois[x,rt]),vrt,length(unique(arois[x,sid])),max(table(arois[x,sid])),length(x))
  })))
  names(currstats)=c(grp,mzmed,"dppm","dppmmin","dppmax",mzmin,mzmax,rt,rtmin,rtmax,"nsid","mndups","n")
  return(currstats)
}
