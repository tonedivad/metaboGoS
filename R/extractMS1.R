### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name refinePIRoi
#' @title refine PI ROI
#' 
#' Refine/integra ROI around spectra  after association of spectra to ROIs from potential precursors
#'
#' @param obj metaboGoS obj
#' @param parDeco deconvolution parameter list
#' @param refFct refining function
#' @param doPlot plotting at each step
#' @param blList list of blank filles xcmsRaw object 
#' @param nSlaves number of slaves for parallel
#' @param useMS1file use the MS1 scan in MS1file
#' @param keepOld keep the former RoiInfos data frame
#' @import foreach
#' @import doParallel
#' @return metaboGoS obj with RoiInfos / PeakInfos /params
#' @export
refinePIRoi<-function(obj,
                       parDeco, #=list(span=7,bw=0.01,minNoiseMS1=10000),
                       refFct=metaboGoS:::.MGrefFct0,
                      doPlot=TRUE,blList=NULL, nSlaves=1,useMS1file=TRUE,keepOld=TRUE){
  
  ## blList=NULL;nSlaves=1;useMS1file=TRUE;doPlot=TRUE
  
  if(!"ddaSet"%in%class(obj)) stop("Error: not a ddaSet object")
  
  strt=as.integer(as.POSIXct( Sys.time() ))
  
  ########  ########  ########
  oparDeco=obj$parDeco
  if(!is.null(oparDeco))  for(i in names(oparDeco)[!names(oparDeco)%in%parDeco]) parDeco[[i]]=oparDeco[[i]]
  
  # blList=list()
  # if(!is.null(blFiles)){
  #   blFiles=blFiles[file.exists(blFiles)]
  #   if(length(blFiles)>0)  for(i in 1:length(blFiles)) blList[[i]]=xcmsRaw(blFiles[i],includeMSn = FALSE)
  # }
  # 
  
  cefi=normalizePath(paste0(obj$File$dirName,"/",ifelse(!is.null(obj$File$MS1file) & useMS1file,obj$File$MS1file,obj$File$fileName)))
  xr=xcmsRaw(cefi,includeMSn = FALSE)
  sc2rt=round(xr@scantime/60,5)
  
  
  cat("Refining/grouping potential PIs from '",unclass(xr@filepath),"'\n",sep="")
  
  psdrt=ifelse(is.null(parDeco$psdrt),NA,parDeco$psdrt)
  if(is.na(psdrt)) psdrt=parDeco$psdrt=quantile(sc2rt/60,.1)

  ### set new pseudo Scan 
  newrtx=(1:ceiling(max(xr@scantime/60/psdrt)))*psdrt
  sc2nrt=apply(abs(outer(newrtx,xr@scantime/60,"-")),2,which.min)
  sc2nrt=cbind(scan=1:length(xr@scantime),nscan=sc2nrt,nrt=round(newrtx[sc2nrt],5))
  ###
  # if(!is.null(blList)) blList=lapply(blList,function(xrb){
  # sc2nrt2=apply(abs(outer(newrtx,xrb@scantime/60,"-")),2,which.min)
  # sc2nrt2=cbind(scan=1:length(xrb@scantime),nscan=sc2nrt2,nrt=round(newrtx[sc2nrt2],5))
  # list(xrb,sc2nrt2)
  # })
  
  minRTwin=ifelse(is.null(parDeco$minRTwin),NA,parDeco$minRTwin)
  if(is.na(minRTwin)) parDeco$minRTwin=minRTwin=quantile(diff(xr@scantime/60),.9)+2.1*psdrt ## set it to near largest delta rt + 2*psdrt
  
  linpardeco=c("ppm","dmz","psdrt","drt","rtlim","minHeightMS1","minNoiseMS1", "minZero","minRTwin","span","bw" )
  lmiss=linpardeco[!linpardeco%in%names(parDeco)]
  if(length(lmiss)) stop(paste0('Missing parameters: ',paste(lmiss,sep=" ")))
  
  ROImat=obj$RoiInfos
  if(keepOld) obj$OldRoiInfos=obj$RoiInfos
  if(is.null(ROImat)) stop('Matrix defing ROI is missing!')
  ########  ########  ########
  
  ###################### Refine frames
  llre=list()
  ll=ROImat$RoiId#[111:120]
  lperc=ll[round(seq(1,length(ll),length=12)[2:11])]
  
  if(nSlaves>1)   nSlaves=max(1, min(nSlaves,detectCores()-1))
  if(nSlaves>1){
    clProc<-makeCluster(nSlaves)
    doParallel::registerDoParallel(clProc)
    cat(" -- registering ",nSlaves," clusters\n",sep="")
  }

   ### Parallele bit
  if(nSlaves>1)
    #llre=foreach(i = ll, .export = fct2exp,.packages = c("igraph","xcms","GRMeta"), .verbose =F)  %dopar%{
      llre=foreach(idx = ll,.packages = c("metaboGoS"), .verbose =FALSE)  %dopar%{
        lmz=range(ROImat[idx,c("mzmin","mzmax")])
        lrt=range(ROImat[idx,c("rtmin","rtmax")])+c(-3,3)*parDeco$psdrt*(parDeco$span+1) ## add 3*span to make use not in the downslope
        eic=GRMeta:::.GRrawMat(xr,mzrange = lmz, rtrange = lrt*60,padsc =TRUE)
        if(sum(!is.na(eic[,"y"]))<2) return(list())
        if(max(eic[,"y"],na.rm=T)<parDeco$minHeightMS1) return(list())
        eicbl=NULL
        if(length(blList)>0){
          eicbl=do.call("rbind",lapply(1:length(blList),function(ibl){
            i=cbind(GRMeta:::.GRrawMat(blList[[ibl]],mzrange = lmz, rtrange = lrt*60+c(-10,10),padsc =T),blid=ibl)
            if(sum(!is.na(i[,"mz"]))<5) return(NULL)
            i
          }))
        }
        eic=cbind(eic,sc2nrt[eic[,"scan"],2:3])
        re=refFct(eic,eicbl,parDeco,idx,doPlot = FALSE) ## 
        if(is.null(re)) return(list())
        re
    }
  ## Serial bit
  if(nSlaves<=1) for(idx in  ll){
    if(idx %in% lperc) cat(idx,"(",which(ll==idx),") ",sep="")
    lmz=range(ROImat[idx,c("mzmin","mzmax")])
    lrt=range(ROImat[idx,c("rtmin","rtmax")])+c(-3,3)*parDeco$psdrt*(parDeco$span+1) ## add 3*span to make use not in the downslope
    eic=GRMeta:::.GRrawMat(xr,mzrange = lmz, rtrange = lrt*60,padsc =T)
    if(sum(!is.na(eic[,"y"]))<2) next
    if(max(eic[,"y"],na.rm=T)<parDeco$minHeightMS1) next
    eicbl=NULL
    if(length(blList)>0){
      eicbl=do.call("rbind",lapply(1:length(blList),function(ibl){
        i=cbind(GRMeta:::.GRrawMat(blList[[ibl]],mzrange = lmz, rtrange = lrt*60+c(-10,10),padsc =T),blid=ibl)
        if(sum(!is.na(i[,"mz"]))<5) return(NULL)
        i
        }))
    }
    eic=cbind(eic,sc2nrt[eic[,"scan"],2:3])
    re=refFct(eic,eicbl,parDeco,idx,doPlot = doPlot) ## change here + update the 
    if(is.null(re)) next
    llre=c(llre,list(re))
  }

  if(nSlaves>1) stopCluster(clProc)

  ## combine results
  llre=llre[sapply(llre,length)>0]
  if(length(llre)==0){
    cat("Something went wrong no acceptable ROI were found!!!")
    return(obj)
  }
  ares=do.call("rbind",llre)
  ares=ares[which(ares$pk.rt>=parDeco$rtlim[2] & ares$pk.rt<=parDeco$rtlim[3]),]
  
  
  ## Peak matrix
  lvar=unique(c("pk.id","cl.id","roi.id",names(ares)[c(grep("pk\\.",names(ares)))]))
  PKmat=ares[,lvar]
  names(PKmat)[1:3]=c("PkId",'ClId','RoiId')
  names(PKmat)=gsub("^pk.","",names(PKmat))
  ldups=names(which(table(PKmat$PkId)>1))
  for(i in ldups) PKmat$PkId[PKmat$PkId==i]=paste0(PKmat$PkId[PKmat$PkId==i],"-",1:sum(PKmat$PkId==i))
  
  ## ROimatrix
  lvar=unique(c("roi.id",names(ares)[grep("^roi",names(ares))]))
  ROImat2=ares[,lvar]
  names(ROImat2)[1]=c('RoiId')
  names(ROImat2)=gsub("^roi.","",names(ROImat2))
  add2roi=do.call("rbind",tapply(1:nrow(PKmat),PKmat$RoiId,
                                 function(x) round(c(max(PKmat$int.ap[x]),max(PKmat$sbr[x]),round(range(PKmat[x,c("mzmin",'mzmax')]),5),
                                 rt=PKmat$rtap[x[which.max(PKmat$int.ap[x])]],range(PKmat[x,c("rtmin",'rtmax')])),4)))
  colnames(add2roi)=c('intensity',"sbr","mzmin","mzmax","rt","rtmin","rtmax")
  ROImat2=cbind(ROImat2[match(rownames(add2roi),ROImat2$RoiId),],add2roi)
  cat(" -- number of Peaks/ROIs after refinement: ",nrow(PKmat),"/",nrow(ROImat2)," [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")
  
  l2k=which(ROImat2[,"intensity"]>=parDeco$minHeightMS1 & abs(ROImat2[,"rtmax"]-ROImat2[,"rtmin"])>=parDeco$minRTwin &
              ROImat2[,"rtmin"]<=max(parDeco$rtlim) & ROImat2[,"rtmax"]>=min(parDeco$rtlim) & ROImat2$coda>0.2 & ROImat2$sbr>=2)
  ROImat2=ROImat2[l2k,,drop=F]
  ROImat2=ROImat2[order(ROImat2[,'mz50'],ROImat2[,"rtmin"]),,drop=F]
  rownames(ROImat2)=ROImat2$RoiId #sprintf("R%.4f@%.1f-%.1f",ROImat2[,"mz50"],ROImat2[,"rtmin"],ROImat2[,"rtmax"])

  PKmat=PKmat[which(PKmat$RoiId%in%ROImat2$RoiId),]
  rownames(PKmat)=PKmat$PkId
 
  ###########################
  ### Associate MS/MS to ROI
  ## !! recompute the exact values of the precursors??
  oMS2ROI=obj$MS2toMS1
  oMS2ROI$OldRoiId=oMS2ROI$RoiId
  llsplit=split(1:nrow(oMS2ROI), ceiling(seq_along(1:nrow(oMS2ROI))/200))
  Roi2sp=do.call("rbind",lapply(llsplit,function(x){
    ddmz=outer(oMS2ROI$mzPrec[x],ROImat2[,"mzmin"],"-")>=(-parDeco$dmz/2) & outer(oMS2ROI$mzPrec[x],ROImat2[,"mzmax"],"-")<=(parDeco$dmz/2) ## slight padding zeros instaeds?
    ddrt=outer(oMS2ROI$rtPrec[x],ROImat2[,"rtmin"],"-")>=-(parDeco$minRTwin*2) & outer(oMS2ROI$rtPrec[x],ROImat2[,"rtmax"],"-")<=(parDeco$minRTwin*2)
    re=which(ddrt & ddmz,arr=T)
    re[,1]=x[re[,1]]
    re}))
  # col: ROImat2, row=oMS2ROI
  nMS2ROI=oMS2ROI[Roi2sp[,"row"],]
  nMS2ROI$RoiId=ROImat2$RoiId[Roi2sp[,"col"]]
  nMS2ROI=nMS2ROI[,which(!names(nMS2ROI)%in%c("RoiPkCl" , "PkId", "OldRoiId"))]
  
  ### Add Peak to spectra
 # oMS2ROI=obj$MS2toMS1=NA
  llsplit=tapply(1:nrow(nMS2ROI),nMS2ROI$RoiId,c)
  Pk2sp=do.call("rbind",lapply(llsplit,function(x){
    iroi=nMS2ROI$RoiId[x][1]
    l=which(PKmat$RoiId==iroi)
    if(length(l)==0) return(NULL)
    re=which(outer(nMS2ROI$rtPrec[x],PKmat$rtmin[l],"-")>=(-parDeco$minRTwin*4) & outer(nMS2ROI$rtPrec[x],PKmat$rtmax[l],"-")<=(parDeco$minRTwin*4),arr=T)
    re[,1]=x[re[,1]]
    re[,2]=l[re[,2]]
    re}))
  # col: peak, row=nMS2ROI
  
  fMS2ROI=nMS2ROI[Pk2sp[,1],c("SpId","RoiId","PInt","yPrec","mzPrec", "rtPrec")]
  fMS2ROI$PkId=PKmat$PkId[Pk2sp[,2]]
  fMS2ROI$ClId=PKmat$ClId[Pk2sp[,2]]
  l2add=which(!(1:nrow(nMS2ROI))%in%Pk2sp[,1])
  if(length(l2add)){
    l2add=nMS2ROI[l2add,c("SpId","RoiId","PInt","yPrec","mzPrec", "rtPrec")]
    l2add$ClId=l2add$PkId=NA
    fMS2ROI=rbind(fMS2ROI,l2add)
  }
  fMS2ROI=fMS2ROI[order(fMS2ROI$mzPrec,fMS2ROI$rtPrec,fMS2ROI$RoiId,fMS2ROI$PkId),]
  rownames(fMS2ROI)=NULL
 
  cat(" --> Final number of ROIs: ",nrow(ROImat2)," assoc. to ",length(unique(fMS2ROI$SpId)),"/",nrow(obj$MS2Infos)," MS/MS\n",
      " --> Final number of peaks: ",nrow(PKmat)," assoc. to ",length(unique(fMS2ROI$SpId[!is.na(fMS2ROI$PkId)])),"/",nrow(obj$MS2Infos)," MS/MS",
      " [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")
  
  
  obj$RoiInfos=ROImat2
  obj$PeakInfos=PKmat
  obj$MS2toMS1=fMS2ROI
  obj$parDeco=parDeco
  class(obj) = unique(append(class(obj), "ddaSet"))
  invisible(obj)
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name extractPIRoi
#' @title Extract ROI around spectra
#' 
#' Fast association of spectra to ROIs from potential precursors
#'
#' @param obj metaboGoS obj
#' @param parDeco deconvolution parameter list
#' @param minPIwindow window size in min. around the spectra
#' @param PIweight weighting
#' @param useMS1file use the MS1 scan in MS1file
###' @param nSlaves number of slaves for parallel
#' @import foreach
#' @import doParallel
#' @return metaboGoS obj with RoiInfos / MS2ROI /params
#' @export
extractPIRoi<-function(obj,
                    parDeco=list(minPurity=c(0.2,0.4),
                         ppm=5,dmz=0.001,psdrt=NA,
                         drt=1,rtlim=c(-Inf,-Inf,Inf,Inf),
                         minHeightMS1=100000,minZero=200),
                    minPIwindow=NA,PIweight=2,useMS1file=TRUE){
  
  ## minPIwindow=0.05;PIweight=2;nSlaves=10
  if(!"ddaSet"%in%class(obj)) stop("Error: not a ddaSet object")

  strt=as.integer(as.POSIXct( Sys.time() ))

  ms2inf=obj$MS2Infos
  cefi=normalizePath(paste0(obj$File$dirName,"/",ifelse(!is.null(obj$File$MS1file) & useMS1file,obj$File$MS1file,obj$File$fileName)))
  xr=xcmsRaw(cefi,includeMSn = FALSE)
  cat("Extract potential PIs from '",unclass(xr@filepath),"' comprising ",nrow(ms2inf)," MS/MS and ",length(xr@scantime)," MS1\n",sep="")
  
  sc2rt=round(xr@scantime/60,5)
  psdrt=ifelse(is.null(parDeco$psdrt),NA,parDeco$psdrt)
  if(is.na(psdrt)) psdrt=parDeco$psdrt=quantile(sc2rt/60,.1)
  
  if(is.na(minPIwindow)) minPIwindow=quantile(diff(xr@scantime/60),.9)+2.1*psdrt ## set it to near largest delta rt + 2*psdrt
  parDeco$minRTwin=minRTwin=minPIwindow
  
  if(any(is.na(ms2inf$WinSize))) stop('Please set the delta m/z around precursor ions')
  
  rtlim=parDeco$rtlim
  rtlim[which(rtlim>max(sc2rt))]=max(sc2rt)
  rtlim[which(rtlim<min(sc2rt))]=min(sc2rt)
  parDeco$rtlim=rtlim
  
  ##### Check
  linpardeco=c("minPurity" ,   "ppm","dmz","psdrt","drt","rtlim","minHeightMS1", "minZero","minRTwin" )
  lmiss=linpardeco[!linpardeco%in%names(parDeco)]
  if(length(lmiss)) stop('Missing parameters:',paste(lmiss,sep=" "))

  ### Get potential precursor within minRTwin of the spectra
  amzrt=list()
  for(ipmz in 1:nrow(ms2inf)){
    linrange=which(abs(sc2rt-ms2inf[ipmz,'RT'])<=minRTwin)
    if(length(linrange)==0) next
    scrange=range(linrange)
    imz=ms2inf[ipmz,]$PrecMZ
    mzrange=imz+c(-1,1)*(ms2inf[ipmz,]$WinSize+parDeco$dmz)
    re=GRMeta:::.GRrawMat(xr,scanrange=scrange,mzrange=mzrange)
    if(nrow(re)==0) next #return(NULL)# return(c(PInt=NA,MZ=NA,SBR=NA,CodaWide=NA))
    we=.MGgetMZweight(re[,"mz"]-imz,PIweight)
    re[,"y"]=re[,"y"]*we
    pint=re[,"y"]/tapply(re[,"y"],re[,"scan"],sum)[as.character(re[,"scan"])]
    if(max(pint)<parDeco$minPurity[1]) next #return(NULL)# return(c(PInt=NA,MZ=NA,SBR=NA,CodaWide=NA))
    re=cbind(i=ipmz,pint=pint,re[,c("mz","rt","y"),drop=F])
    amzrt[[ipmz]]=re[pint>=parDeco$minPurity[1],,drop=F]
  }
  amzrt=do.call("rbind",amzrt)
  cat(" -- ",nrow(amzrt)," ions found: ",length(unique(amzrt[,1]))," spectra out of ",nrow(ms2inf),"\n Purity summary:\n",sep="")
  print(summary(amzrt[,2]))
  
  ### merge them
  llppm<-GRMeta:::.GRsplistMZ(amzrt[,"mz"],dppm = parDeco$ppm,dmz=parDeco$dmz,typ="min") ## min distance b/w dppm/dmz
  ROImat=do.call("rbind",lapply(llppm,function(x) c(median(amzrt[x,"mz"]),range(amzrt[x,"mz"]),round(range(amzrt[x,"rt"]),4))))
  colnames(ROImat)=c("mz","mzmin","mzmax","rtmin","rtmax")
  ROImat=ROImat[order(ROImat[,1]),,drop=F]
  cat(" -- reduced to ",nrow(ROImat)," ROIs [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")
  
  
  ### merge wide m/z/RT windows
  ROImat1=lapply(1:nrow(ROImat),function(ix){
    xroi=ROImat[ix,]
    lmz=range(c(xroi[c( "mz","mz","mzmin", "mzmax")]*(1+1*c(-2,2,-1,1)*parDeco$ppm*10^-6),xroi[c("mzmin", "mzmax")]+c(-1,1)*parDeco$dmz))+c(-1.1,1.1)*parDeco$dmz ## very large window!!
    lrt=range(xroi[c( "rtmin","rtmax")])+c(-1,1)*parDeco$drt
    ieic=xcms::rawEIC(xr, mzrange =lmz, rtrange = lrt*60)
    l=which(ieic$intensity>=parDeco$minZero)
    ll=GRMeta:::.GRsplist(sc2rt[l],l,d=parDeco$drt*2+2.1*psdrt)
    do.call("rbind",lapply(ll,function(x){ 
      newx=sc2rt[range(ieic$scan[x])]+c(-1,1)*(parDeco$drt)
      #x2=x[ieic$scan[x]>=rangeSc[2] & ieic$scan[x]<=rangeSc[3]]
      #if(length(x2)<1) return(c(ieicoi[1],lmz[1],lmz[2],newx,0,0,0))
      aprt=ieic$scan[x[which.max(ieic$intensity[x])]]
      return(c(xroi[1],lmz[1],lmz[2],newx,diff(newx),max(ieic$intensity[x]),sc2rt[aprt]))
    }))
  })
  ROImat1=do.call("rbind",ROImat1)
  colnames(ROImat1)=c(colnames(ROImat),"drt","intensity","rt")
  
  l2k=which(ROImat1[,"intensity"]>=parDeco$minHeightMS1 & ROImat1[,"drt"]>=minRTwin & ROImat1[,"rtmin"]<=rtlim[4] & ROImat1[,"rtmax"] >= rtlim[1])
  ROImat1=ROImat1[l2k,]
  ROImat1[ROImat1[,"rtmin"]<=rtlim[1],"rtmin"]=rtlim[1]
  ROImat1[ROImat1[,"rtmax"]>=rtlim[4],"rtmax"]=rtlim[4]
  
  ### Combine dupl -> remerge retention time -> could be improved!?!
  llover=GRMeta:::.GRisover(ROImat1[,"mzmin"],ROImat1[,"mzmax"],retOne = T)
  ROImat2=do.call("rbind",lapply(llover,function(x){
    if(length(x)==0) return(ROImat1[x,,drop=F])
    c(mean(ROImat1[x,"mz"]),min(ROImat1[x,"mzmin"]),max(ROImat1[x,"mzmax"]),min(ROImat1[x,"rtmin"]),max(ROImat1[x,"rtmax"]),0,
      max(ROImat1[x,"intensity"]),ROImat1[x[which.max(ROImat1[x,"intensity"])],"rt"])
  }))
  colnames(ROImat2)=colnames(ROImat1)
  ROImat2[,"drt"]=ROImat2[,"rtmax"]-ROImat2[,"rtmin"]
  l2k=which(ROImat2[,"intensity"]>=parDeco$minHeightMS1 & ROImat2[,"drt"]>=minRTwin & ROImat2[,"rtmin"]<=rtlim[4] & ROImat2[,"rtmax"]>=rtlim[1])
  ROImat2=ROImat2[l2k,]
  
  cat(" -- found: ",nrow(ROImat2)," potential ROIs [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")

    ### Associate MS/MS to ROI
  llsplit=split(1:nrow(amzrt), ceiling(seq_along(1:nrow(amzrt))/200))
  Roi2sp=do.call("rbind",lapply(llsplit,function(x){
    ddmz=outer(amzrt[x,"mz"],ROImat2[,"mzmin"],"-")>=(-parDeco$dmz/2) & outer(amzrt[x,"mz"],ROImat2[,"mzmax"],"-")<=parDeco$dmz/2
    ddrt=outer(amzrt[x,"rt"],ROImat2[,"rtmin"],"-")>=0 & outer(amzrt[x,"rt"],ROImat2[,"rtmax"],"-")<=0
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
  ROImat3=data.frame(RoiId=NA,ROImat2[lurois,,drop=FALSE])
  ROImat3$RoiId=rownames(ROImat3)=sprintf("R%.4f@%.1f-%.1f",ROImat3[,"mz"],ROImat3[,"rtmin"],ROImat3[,"rtmax"])
  
  SpId2ROI=data.frame(SpId=ms2inf$SpId[Roi2sp[,"i"]],
                      RoiId=ROImat3$RoiId[Roi2sp[,"col"]],
                     PInt=round(Roi2sp[,"pint"],4),
                     yPrec=round(Roi2sp[,"y"],6),
                     mzPrec=round(Roi2sp[,"mz"],6),
                     rtPrec=round(Roi2sp[,"rt"],6))

  cat(" --> Final number of ROIs: ",nrow(ROImat3)," assoc. to ",length(unique(SpId2ROI$SpId)),"/",nrow(ms2inf),
      " MS/MS [+",as.integer(as.POSIXct( Sys.time() ))-strt,"sec.]\n",sep="")

  obj$RoiInfos=ROImat3
  obj$MS2toMS1=SpId2ROI
  obj$parDeco=parDeco
  class(obj) = unique(append(class(obj), "ddaSet"))
  invisible(obj)
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
    colnames(winmat)=c( "rtmin","rtmax","mz50","mz10","mz90","intensity" ,"mzap","mzmin","mzmax" )
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
.MGgetConsecutive<-function(xeic,ncons=3,ppm=5,dmz=0.001,mindmz=dmz,minRTwin=0.01,minNoise=10000,minHeight=minNoise*2){
#  re=.MGgetConsecutive(xeic,ncons=2,ppm=parDeco$ppm,dmz=parDeco$dmz,mindmz,minRTwin=minRTwin,minNoise=parDeco$minNoiseMS1,minHeight=0)
# ncons=2;ppm=parDeco$ppm;dmz=parDeco$dmz;mindmz;minRTwin=minRTwin;minNoise=parDeco$minNoiseMS1;minHeight=0
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

