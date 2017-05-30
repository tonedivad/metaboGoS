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
#' @param noise.local y/baseline
#' @param span bandwidth
#' @param snr.thresh drt
#' @param minNoise min noise for integrating
#' @param v2alim merge close apices if valley to apex less than 0.9
#' @return a brand new eic
#' @export
.MGsimpleIntegr<-function (x, y, noise.local, span = 5, snr.thresh=NA,minNoise=min(y,noise.local)*1.01,v2alim=0.8,span2=span+2){
  index <- GRMeta:::.GRmsExtrema(y, span = span)
  nvar <- length(x)
  index.min = index$index.min
  index.max = index$index.max
  imax <- which(index.max)
  if(!is.na(snr.thresh)) good.snr = (y[imax]/noise.local[imax])>=snr.thresh else good.snr = (noise.local[imax])
  tick.loc <- imax[good.snr]
  if (length(tick.loc) == 0)    return(data.frame())
  tick.left <- tick.right <- rep(1, length(tick.loc))
  for (i in 1:length(tick.loc)) {
    tick.left[i] <- max(which(index.min[1:tick.loc[i]]))
    tick.right[i] <- min(which(index.min[tick.loc[i]:length(x)])) + tick.loc[i] - 1
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
  
  ### 
  l2m=NULL
  if(!is.na(v2alim)){
    releft=y[pks[,2]]/y[pks[,1]]
    reright=y[pks[,3]]/y[pks[,1]]
    releft[1]=reright[length(reright)]=0
    l2mle=which(releft>=v2alim & abs(pks[,2]-pks[,1])<=(span2))
    if(length(l2mle)>0) l2m=cbind(l2mle-1,l2mle)
    l2mri=which(reright>=v2alim & abs(pks[,3]-pks[,1])<=(span2))
    if(length(l2mri)>0) l2m=rbind(l2m,cbind(l2mri,l2mri+1))
  }
  ##
  if(!is.null(l2m)){
    l2m=.GRmergellx(lapply(1:nrow(l2m),function(x) l2m[x,]))
    pks=rbind(pks,do.call("rbind",lapply(l2m,function(ix) c(pks[ix,1][which.max(y[pks[ix,1]])],range(pks[ix,2:3])))))
    pks=pks[-unlist(l2m),,drop=F]
    pks=pks[order(pks[,1]),,drop=F]
  }
  
  
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
                   pk.snr = noise.local[pks[, 1]],
                   pk.int = y[pks[, 1]], 
                   pk.span = x[pks[, 3]] - x[pks[, 2]])
  pks$pk.cl = lcl
  
  
  return(pks)
}

######  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ######## 
#'
#' rewritting of ksmooth to fill-in missing scan/uneually spaced rt
#'
#' @param x retention time
#' @param y intensity
#' @param noise.local +1/0/-1
#' @param span bandwidth
#' @param span2 merging span
#' @param snr.thresh drt
#' @param minNoise min noise for integrating
#' @param v2alim merge close apices if valley to apex less than 0.9
#' @return a brand new eic
#' @export
.MGsimpleIntegr2<-function (x, y, bsl,bslscore,snr.thresh=2, span = 5,minNoise=min(y,bsl)*1.01,v2alim=0.8,span2=2*span+1-2){
  
  index <- GRMeta:::.GRmsExtrema(y, span = span)
  nvar <- length(x)
  index.min = index$index.min
  index.max = index$index.max
  imax <- which(index.max)
  tick.loc <-imax[((y/bsl)>=snr.thresh & bslscore>1)[imax]]
  tick.loc2 <-imax[(bslscore[imax]>1)]
  tick.loc2=tick.loc2[!tick.loc2%in%tick.loc]
  if (length(tick.loc) == 0) return(data.frame())
  
  ### remove valleys to close to each other 
  ival<- which(index.min)
  l2merge=lapply(2:length(ival),function(i){
  #  if(!any(index.min[tick.loc[i-1]:tick.loc[i]])) return(c(tick.loc[i-1],tick.loc[i]))
    if(max(y[c(ival[i],ival[i-1])]/max(y[ival[i]:ival[i-1]]))>=v2alim & (ival[i]-ival[i-1])<span+3) return(c(ival[i-1],ival[i]))
    c()})
  l2merge=l2merge[sapply(l2merge,length)>1]
  if(length(l2merge)>0){
    l2merge=.GRmergellx(l2merge)
    for(i in l2merge){
      val2rm=i[-which.min(bslscore[i])]
  #    pk2rm=which(tick.loc>=min(val))
      index.min[val2rm]=FALSE
    }  
  }
  
  ### check if tick.lock have index.min
  if(length(tick.loc)>1){
    l2merge=lapply(2:length(tick.loc),function(i){
      if(!any(index.min[tick.loc[i-1]:tick.loc[i]])) return(c(tick.loc[i-1],tick.loc[i]))
      if(max(min(y[tick.loc[i]:tick.loc[i-1]])/y[c(tick.loc[i],tick.loc[i-1])])>=v2alim & (tick.loc[i]-tick.loc[i-1])<2*span) return(c(tick.loc[i-1],tick.loc[i]))
      c()})
    l2merge=l2merge[sapply(l2merge,length)>1]
    if(length(l2merge)>0){
      l2merge=.GRmergellx(l2merge)
      for(i in l2merge){
        tick.loc=tick.loc[!tick.loc%in%i[-which.max(y[i])]]
        index.min[min(i):max(i)]=FALSE
      }
    }
  }
  
  ## if no valleys outside the tick.loc
  if(min(tick.loc)<min(which(index.min))) index.min[max(which(rev(diff(y[rev(1:min(tick.loc))]))>0)-1,1,na.rm=T)]=TRUE
  if(max(tick.loc)>max(which(index.min))) index.min[max(max(tick.loc)+which(diff(y[max(tick.loc):length(y)])>0)[1],length(y),na.rm=T)]=TRUE
  
  ## if no valleys between tick.loc??
  
  
  # 
  # lmins=unique(c(1,which(index.min & bslscore< -1),length(y)))
  # lmins=lapply(2:length(lmins),function(i) lmins[i-1]:lmins[i])
  # lmins=lmins[sapply(lmins,function(i) any(i%in%tick.loc))]
  # 
  # 
  tick.left <- tick.right <- rep(NA, length(tick.loc))
  ## check left side
  for (i in order(-bslscore[tick.loc],-y[tick.loc])) {
    llv=rev(which(index.min[1:tick.loc[i]]))
    if(i>1) llv=llv[llv>max(tick.loc[1:(i-1)])]
    if(length(llv)==1){tick.left[i]=llv;next}
    if(any(bslscore[llv]< -1)) llv=llv[1:which(bslscore[llv]< -1)[1]] ## pick the closest well below bsl
    if(length(llv)==1){tick.left[i]=llv;next}
    
    while(any(abs(diff(llv))<(1.5*span+1))) llv=llv[-which.min(abs(diff(llv)))]
    if(length(llv)==1){tick.left[i]=llv;next}
    # while(any(abs(diff(llv))>(3*span+1))) llv=llv[-(which.max(abs(diff(llv)))+1)]
    # if(length(llv)==1){tick.left[i]=llv;next}
    while(any(diff(y[llv])>0))  llv=llv[1:which(diff(y[llv])>0)[1]]
    if(length(llv)==1){tick.left[i]=llv;next}
    
    if(any(bslscore[llv]< 0)) llv=llv[1:which(bslscore[llv]< 0)[1]] 
    if(length(llv)==1){tick.left[i]=llv;next}
    if(any(bslscore[llv]< 1)) llv=llv[1:which(bslscore[llv]< 1)[1]] ## pick the closest well below bsl
    tick.left[i]=llv[which.min(bslscore[llv])]
    
  }
  
  ## check right side
  for (i in order(-bslscore[tick.loc],-y[tick.loc])) {
    llv=which(index.min[tick.loc[i]:length(x)])+ tick.loc[i] - 1
    if(i<length(tick.loc)) llv=llv[llv<min(tick.loc[(i+1):length(tick.loc)])]
    if(length(llv)==1){tick.right[i]=llv;next}
    if(any(bslscore[llv]< -1)) llv=llv[1:which(bslscore[llv]< -1)[1]] ## pick the closest well below bsl
    if(length(llv)==1){tick.right[i]=llv;next}
    while(any(diff(llv)<(1.5*span+1))) llv=llv[-which.min(diff(llv))]
    if(length(llv)==1){tick.right[i]=llv;next}
    # while(any(diff(llv)>(3*span+1))) llv=llv[-(which.max(diff(llv))+1)]
    # if(length(llv)==1){tick.right[i]=llv;next} 
    while(any(diff(y[llv])>0))  llv=llv[1:which(diff(y[llv])>0)[1]] ## valley must go down
    if(length(llv)==1){tick.right[i]=llv;next}
    if(any(bslscore[llv]< 0)) llv=llv[1:which(bslscore[llv]< 0)[1]]
    if(length(llv)==1){tick.right[i]=llv;next}
    if(any(bslscore[llv]< 1)) llv=llv[1:which(bslscore[llv]< 1)[1]]
    if(length(llv)==1){tick.right[i]=llv;next}
    tick.right[i]=llv[which.min(bslscore[llv])]
    
  }
  
  # 
  #   llright=
  #   if(i<length(tick.loc)) llright=llright[llright<min(tick.loc[-(1:i)])]
  #    ## pick the closest well below bsl
  #   # tick.left[i] <- max()
  #   # tick.right[i] <- min() 
  # }
  # 
  ## fix left and right
  if (tick.right[length(tick.loc)] == 1) 
    tick.right[length(tick.loc)] <- length(x)
  # if (tick.rightn[length(tick.loc)] == 1) 
  #   tick.rightn[length(tick.loc)] <- length(x)
  if (tick.left[1] == 1) 
    tick.left[1] <- which.min(y[1:tick.loc[1]] != minNoise)
  # if (tick.leftn[1] == 1) 
  #   tick.leftn[1] <- which.min(y[1:tick.loc[1]] != minNoise)
  
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
  
  ### 
  l2m=NULL
  if(!is.na(v2alim)){
    releft=y[pks[,2]]/y[pks[,1]]
    reright=y[pks[,3]]/y[pks[,1]]
    releft[1]=reright[length(reright)]=0
    l2mle=which(releft>=v2alim & abs(pks[,2]-pks[,1])<=(span2))
    if(length(l2mle)>0) l2m=cbind(l2mle-1,l2mle)
    l2mri=which(reright>=v2alim & abs(pks[,3]-pks[,1])<=(span2))
    if(length(l2mri)>0) l2m=rbind(l2m,cbind(l2mri,l2mri+1))
  }
  ##
  if(!is.null(l2m)){
    l2m=.GRmergellx(lapply(1:nrow(l2m),function(x) l2m[x,]))
    pks=rbind(pks,do.call("rbind",lapply(l2m,function(ix) c(pks[ix,1][which.max(y[pks[ix,1]])],range(pks[ix,2:3])))))
    pks=pks[-unlist(l2m),,drop=F]
    pks=pks[order(pks[,1]),,drop=F]
  }
  
  
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
                   pk.snr = (y/bsl)[pks[, 1]],
                   pk.int = y[pks[, 1]], 
                   pk.span = x[pks[, 3]] - x[pks[, 2]])
  pks$pk.cl = lcl
  
  
  return(pks)
}

#######  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ######## 
#'
#' naive integration eic with fill-in missing scan/scna -> first screening of PI must be rewordered to fit last integration step
#'
#' @param xeic EIC
#' @param drt drt
#' @param span bandwidth
#' @param minRTwin baseline
#' @param minNoise min noise for integrating
#' @param minHeight min noise for integrating
#' @param sbslr ratio sample to baseline
#' @return a brand new eic
#' @export
.MGintegrateEIC<-function(xeic,drt,span=5,bw=drt*span*2,minNoise=10000,minHeight=NA,sbslr=2){
  
  # drt = parDeco$psdrt;span=parDeco$span;bw=parDeco$bw;minNoise=parDeco$minNoiseMS1;minHeight=parDeco$minHeightMS1
  nspan=floor(span/2)*2+1
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
  y=newx$y
  y=c(rep(y[1],5*nspan+1),y,rep(rev(y)[1],5*nspan+1))
  
  bsl=GRMeta:::.GRbslrf(1:length(y),y,NoXP = NULL)
  bsl$fit[bsl$fit<minNoise]=minNoise
  bslscore <- (y - bsl$fit)/max(bsl$sigma, 10^-3)
  bslscore[which(abs(bslscore) > 10)] = sign(bslscore[which(abs(bslscore) > 10)]) * 10
  newx$bslc=bslscore[nspan*5+1+(1:length(newx$y))]
  newx$bsl=bsl$fit[nspan*5+1+(1:length(newx$y))]
  
  # x=newx$x;y=newx$y;noise.local =bsl$fit;snr.thresh = 2;span=11
  
  pks=.MGsimpleIntegr(newx$x,newx$y,noise.local =newx$bsl,snr.thresh = sbslr,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8)
  if(nrow(pks)==0) return(list(NULL,newx))
  ## reduced pks
  lineic=data.frame(do.call("rbind",lapply(1:nrow(pks),function(y){
    l=which(xeic[,"rt"]>=pks$pk.left[y] & xeic[,"rt"]<=pks$pk.right[y])
    c(xeic[l[which.max(xeic[l,"y"])],c("mz","rt","y")],range(xeic[l,"mz"],na.rm=T),
      Area=unname(.GRgetArea(xeic[l,"rt"],xeic[l,"y2"]-10000)[1]),scan=sum(!is.na(xeic[l,"mz"])))
  })))
  names(lineic)=c("ap.mz","ap.rt","ap.int","pk.mzmin","pk.mzmax","pk.area","pk.nsc")
  
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
