#######  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ######## 
#'
#' rewritting of ksmooth to fill-in missing scan/uneually spaced rt
#'
#' @param x retention time
#' @param y intensity
#' @param missc window size in min. around the spectra
#' @param bw bandwidth
#' @param drt drt
#' @return a brand new eic
#' @export
.MGdoksmooth<-function (x,y,missc, bw,drt=NA) 
{
  y2 = (missc!=0) * 1
  lmiss = which(missc)
  if (length(lmiss) > 0) {
    y2 = c(y[1], y, y[length(y)])
    rmse = 10^12
    y2[lmiss + 1] = (y2[lmiss] + y2[lmiss + 2])/2
    y2 = y2[-c(1, length(y2))]
    del = 10
    it = 0
    while (del > 10^-4) {
      it = it + 1
      kres = ksmooth(x, y2, kernel = "normal", bandwidth = bw, x.p = x)
      del = abs(rmse - mean((y2[lmiss] - kres$y[lmiss])^2))
      rmse = mean((y2[lmiss] - kres$y[lmiss])^2)
      y2[lmiss] = kres$y[lmiss]
    }
    y = y2
  }
  kres = ksmooth(x, y, kernel = "normal", bandwidth = bw)
  
  if(!is.na(drt)) kres=approx(kres$x,kres$y,seq(min(x),max(x),drt))
  
  return(kres)
}

#######  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ######## 
#'
#' rewritting of ksmooth to fill-in missing scan/uneually spaced rt
#'
#' @param x retention time
#' @param y intensity
#' @param noise.local baseline
#' @param span bandwidth
#' @param snr.thresh drt
#' @param minNoise min noise for integrating
#' @return a brand new eic
#' @export
.MGsimpleIntegr<-function (x, y, noise.local, span = 5, snr.thresh = 2,minNoise=min(y,noise.local)*1.01){
  index <- GRMeta:::.GRmsExtrema(y, span = span)
  nvar <- length(x)
  index.min = index$index.min
  index.max = index$index.max
  imax <- which(index.max)
  snr = y[imax]/noise.local[imax]
  good.snr <- (snr >= snr.thresh)
  snr <- snr[good.snr]
  tick.loc <- imax[good.snr]
  if ((npeak = length(tick.loc)) == 0) 
    return(data.frame())
  tick.left <- tick.right <- rep(1, length(tick.loc))
  for (i in 1:length(tick.loc)) {
    tick.left[i] <- max(which(index.min[1:tick.loc[i]]))
    tick.right[i] <- min(which(index.min[tick.loc[i]:length(x)])) + 
      tick.loc[i] - 1
  }
  if (tick.right[length(tick.loc)] == 1) 
    tick.right[length(tick.loc)] <- length(x)
  if (tick.left[1] == 1) 
    tick.left[1] <- which.min(y[1:tick.loc[1]] != minNoise)
  
  for (ix in 1:length(tick.loc)) {
    l = tick.left[ix]:tick.right[ix]
    tick.left[ix] = max(tick.left[ix], min(l[which(y[l] >minNoise)]) - 1)
    tick.right[ix] = min(tick.right[ix], max(l[which(y[l] >minNoise)]) + 1)
  }
  lpks = lapply(unique(tick.left), function(x) which(tick.left == x))
  
  pks = do.call("rbind", lapply(lpks, function(ix) {
    if (length(ix) > 1) 
      return(c(tick.loc = tick.loc[ix[which.max(y[tick.loc[ix]])]], 
               tick.left = min(tick.left[ix]), tick.right = max(tick.right[ix])))
    c(tick.loc = tick.loc[ix], tick.left = tick.left[ix], tick.right = tick.right[ix])
  }))
  
  ### comp cluster of peaks
  icl = 0
  lcl = l = c()
  for (ix in 1:nrow(pks)) {
    if (any((pks[ix, 3]:pks[ix, 2]) %in% l)) {
      l = c(l, pks[ix, 3]:pks[ix, 2])
    }
    else {
      l = pks[ix, 3]:pks[ix, 2]
      icl = icl + 1
    }
    lcl = c(lcl, icl)
  }
  pks = data.frame(pk.loc = x[pks[,1]],
                   pk.left = x[pks[, 2]],
                   pk.right = x[pks[,3]],
                   pk.snr = y[pks[, 1]]/noise.local[pks[, 1]],
                   pk.int = y[pks[, 1]], 
                   pk.span = x[pks[, 3]] - x[pks[, 2]])
  pks$pk.cl = lcl
  
  
  return(pks)
}
#######  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ######## 
#'
#' naive integration eic with fill-in missing scan/scna 
#'
#' @param xeic EIC
#' @param drt drt
#' @param span bandwidth
#' @param minRTwin baseline
#' @param minNoise min noise for integrating
#' @param minHeight min noise for integrating
#' @return a brand new eic
#' @export
.MGintegrateEIC<-function(xeic,drt,span=5,bw=drt*span*2,minNoise=10000,minHeight=NA){
  
  # drt = parDeco$psdrt;span=parDeco$span;bw=parDeco$bw;minNoise=parDeco$minNoiseMS1;minHeight=parDeco$minHeightMS1
  
  lnnas=which(!is.na(xeic[,"y"]) &  xeic[,"y"]>=minNoise)
  segpks=GRMeta:::.GRsplist(xeic[lnnas,"rt"],lnnas,d=span*drt/1.99,ismass = F) ## allows 2 missing scan if delta scan large due to a lot of 
  xrt=t(sapply(segpks,function(y) range(xeic[y,"rt"])))
  xrt[,1]=xrt[,1]-bw/2
  xrt[,2]=xrt[,2]+bw/2
  if(nrow(xrt)>1){
    ll2merge=GRMeta:::.GRisover(xrt[,1],xrt[,2],T,bw-drt)
    xrt=do.call("rbind",lapply(ll2merge,function(x) range(xrt[x,])))
  }
  xeic=cbind(xeic,"isin"=0)
  for(i in 1:nrow(xrt)) xeic[xeic[,"rt"]>=xrt[i,1] & xeic[,"rt"]<=xrt[i,2],"isin"]=1
  xeic=cbind(xeic,y2=xeic[,"y"])
  xeic[which(xeic[,"isin"]==0 | xeic[,"y"]<=minNoise),"y2"]=minNoise
  if(is.na(xeic[1,"y2"])) xeic[1,"y2"]=minNoise
  if(is.na(xeic[nrow(xeic),"y2"])) xeic[nrow(xeic),"y2"]=minNoise
  
  ## Pseudo filling
  xeic=cbind(xeic,toimp=0)
  lnna=which(is.na(xeic[,"y2"]))
  if(length(lnna)){
    l2rep=GRMeta:::.GRsplist(lnna,lnna,d=1.1)
    for(y in l2rep){
      df=data.frame(Y=log(xeic[range(y)+c(-1,1),"y2"]),RT=xeic[range(y)+c(-1,1),"rt"])
      xeic[y,"y2"]=exp(predict(lm(Y~RT,df),newdata=data.frame(RT=xeic[y,"rt"])))
      xeic[y,"toimp"]=1
    }
  }
  
  newx=.MGdoksmooth(x=xeic[,"rt"],y=xeic[,"y2"],missc = xeic[,"toimp"]==1,bw=bw,drt=drt)
  bsl=GRMeta:::.GRbslrf(newx$x,newx$y,NoXP = NULL)
  bsl$fit[bsl$fit<minNoise]=minNoise
  newx$bsl=bsl$fit
  bslscore <- (newx$y - newx$bsl)/max(bsl$sigma, 10^-3)
  bslscore[which(abs(bslscore) > 10)] = sign(bslscore[which(abs(bslscore) > 10)]) * 10
  newx$bslc=bslscore
  
  # x=newx$x;y=newx$y;noise.local =bsl$fit;snr.thresh = 2;span=11
  
  pks=.MGsimpleIntegr(newx$x,newx$y,noise.local =newx$bsl,snr.thresh = 2,span=floor(span/2)*2+1,minNoise = minNoise*1.01)
  if(nrow(pks)==0) return(list(NULL,newx))
  ## reduced pks
  lineic=data.frame(do.call("rbind",lapply(1:nrow(pks),function(y){
    l=which(xeic[,"rt"]>=pks$pk.left[y] & xeic[,"rt"]<=pks$pk.right[y])
    c(xeic[l[which.max(xeic[l,"y"])],c("mz","rt","y")],Area=unname(.GRgetArea(xeic[l,"rt"],xeic[l,"y2"]-10000)[1]))
  })))
  names(lineic)=c("ap.mz","ap.rt","ap.int","pk.area")
  
  lineic$pk.areasm=unname(sapply(1:nrow(pks),function(y){
    l=which(newx$x>=pks$pk.left[y] & newx$x<=pks$pk.right[y])
    .GRgetArea(newx$x[l],newx$y[l])[1]-.GRgetArea(newx$x[l],newx$bsl[l])[1]
  }))
  
  pks=cbind(pks,lineic)
  if(is.na(minHeight)) return(list(pks,newx))
  l2k=which(tapply(pks$ap.int,pks$pk.cl,max)>=minHeight)
  pks=pks[pks$pk.cl%in%l2k,,drop=F]
  if(nrow(pks)==0) return(list(NULL,newx))
  pks$pk.cl=as.numeric(factor(pks$pk.cl))
  return(list(pks,newx))
}



#######  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ######## 
