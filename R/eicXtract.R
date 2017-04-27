### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name eicXtract
#' @title refine PI ROI
#' 
#' Refine/integra ROI around spectra  after association of spectra to ROIs from potential precursors
#'
#' @param infile metaboGoS obj
#' @param eicmat eicmat
#' @param addppm add ppm 
#' @param refFct refining function
#' @param drt plotting at each step
#' @param nSlaves number of slaves for parallel
#' @param save numbl
#' @param outfile n
#' @import foreach
#' @import doParallel
#' @return metaboGoS obj with RoiInfos / PeakInfos /params
#' @export
eicXtract<-function(infile,eicmat,addppm=NA,refFct=metaboGoS:::.MGfill,drt,nSlaves=1,save=FALSE,outfile=NULL,...){
  
print(infile)
#drt=0.00446
xr=xcmsRaw(infile)
if(!is.na(addppm)) xr=.GRcorrRawPPM(xr,addppm)
newrtx=(1:ceiling(max(xr@scantime/60/drt)))*drt
sc2nrt=apply(abs(outer(newrtx,xr@scantime/60,"-")),2,which.min)
sc2nrt=cbind(scan=1:length(xr@scantime),nscan=sc2nrt,nrt=round(newrtx[sc2nrt],5))
#sc2nrt=round(newrtx[sc2nrt],6)

allxeic=list()
ll=eicmat$RoiId#[101:150]
lperc=ll[round(seq(1,length(ll),length=12)[2:11])]

if(nSlaves>1)   nSlaves=max(1, min(nSlaves,detectCores()-1))
if(nSlaves>1){
  clProc<-makeCluster(nSlaves)
  doParallel:::registerDoParallel(clProc)
  cat(" -- registering ",nSlaves," clusters\n",sep="")
}

### Parallele bit
if(nSlaves>1){
  #llre=foreach(i = ll, .export = fct2exp,.packages = c("igraph","xcms","GRMeta"), .verbose =F)  %dopar%{
  allxeic=foreach(iroi = ll,.packages = c("metaboGoS"), .verbose =F)  %dopar%{
    i=which(eicmat$RoiId==iroi)
    xeic=.GRrawMat(xr,mzrange=range(eicmat[i,c("mzmin","mzmax")]),rtrange = (range(eicmat[i,c("rtmin","rtmax")])+c(-2,2)*drt)*60,padsc = T) # extrat pad for the newscan
    if(all(is.na(xeic[,"y"]))) return(NULL)
    xeic=cbind(xeic,sc2nrt[xeic[,"scan"],2:3])
    df=data.frame(do.call("cbind",.MGfill(xeic,drt=drt,...)))
    return(list(iroi,df))
  }
  allxeic=allxeic[sapply(allxeic,length)==2]
  nallxeic=sapply(allxeic,function(x) x[[1]]) 
  allxeic=lapply(allxeic,function(x) x[[2]])
  names(allxeic)=nallxeic
  
}
## Serial bit
if(nSlaves<=1) for(iroi in ll){
  if(iroi%in%lperc) cat(iroi," ")
  i=which(eicmat$RoiId==iroi)
  rtrange=(range(eicmat[i,c("rtmin","rtmax")])+c(-3,3)*drt)
xeic=.GRrawMat(xr,mzrange=range(eicmat[i,c("mzmin","mzmax")]),rtrange = rtrange*60,padsc = T) # extrat pad for the newscan
if(all(is.na(xeic[,"y"]))) next
xeic=cbind(xeic,sc2nrt[xeic[,"scan"],2:3])
#if(any(table(xeic[,"scan"])>1)) print(c(sum(table(xeic[,"scan"])>1),i))
df=data.frame(do.call("cbind",.MGfill(xeic,drt=drt,...)))
#df=data.frame(do.call("cbind",.MGfill(xeic,drt=drt,span=7)))
allxeic[[iroi]]=df
}

if(nSlaves>1) stopCluster(clProc)
if(save & is.null(outfile)) save(file=paste0(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", infile),".rda"),allxeic,sc2nrt)
if(!is.null(outfile)) save(file=outfile,allxeic,sc2nrt)
invisible(allxeic)

}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' 
#' ..
#' 
#' @param xeic EIC
#' @param drt weight
#' @param span weight
#' @param bw weight
#' @param minNoise weight
#' @keywords internal
#' 
#' @export
.MGfill<-function(xeic, drt, span = 5, bw = drt * span * 2, minNoise = 5000){
  
  # span = 5;bw = drt * span * 2;minNoise = 10000;minHeight = NA
  
  
  lnnas = which(!is.na(xeic[, "y"]) & xeic[, "y"] >= minNoise)
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

