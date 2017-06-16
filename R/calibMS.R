#' @name calibMS1
#' 
#' @title Convert CFM file to MGF
#' some stufff
#'
#' @param filein Path to the input file or xcmsRaw object
#' @param lmzcalib Calibration masses or method name in BGIons
#' @param dppm Mass error
#' @param bw Bandwidth to group masses
#' @param minpts Minimum of points
#' @param bkpoint breakpoint mass to calculate robust regression
#' @param calInt calibration interval to compute 
#' @param doPlot should plot
#' @import segmented
#' @return list of matches and calibration function/paramters
#' @export
calibMS1<-function(filein,lmzcalib,dppm=20,bw=2.5,minpts=31,bkpoint=median(calInt),calInt=c(70,500),typ=1,doPlot=T,retMat=F){
  
  ## dppm=20;bw=2.5;minpts=31;bkpoint=median(calInt);calInt=c(70,500);typ=1;doPlot=T;retMat=F
  
  ###
  .infct<-function(matdf,calpars,calfct,typ=3){
    fail3=FALSE
    if(typ==3){
      fct2<-function(x,a,b,brk)  b*(x >= brk) + (b - (a/2) * (brk-x)^2) * (x < brk)
      y=matdf$dppm[matdf$out]
      x=matdf$calmz[matdf$out]
      b=median(y[x>calpars["brk"]])
      x=x-(calpars["brk"]*1.2) ## a few more points above the estimatated breakpoint
      y=y-b
      a=abs(coef(lm(y[x<0]~I(x[x<0]^2)))[2])*sign(mean(y[x<0]))*-2 ## make sure this is the correct sign
      m <- try(nls(dppm ~ fct2(calmz,a,b,brk),
                   start = list(a = a, b = b,brk=calpars["brk"]),
                   data=matdf,subset = out),TRUE)
      if(!"try-error"%in%class(m)){
        calpars<-c(coef(m),low=calInt[1],high=calInt[2])
        calfct<-function(x,pars,round=3){
          pr= pars["b"]*(x >=  pars["brk"]) + (pars["b"] - (pars["a"]/2) * ( pars["brk"]-x)^2) * (x <  pars["brk"])
          x2=pars[c("low","high")]
          pr2= pars["b"]*(x2 >=  pars["brk"]) + (pars["b"] - (pars["a"]/2) * ( pars["brk"]-x2)^2) * (x2 <  pars["brk"])
          pr[x<pars["low"]]=pr2[1]
          pr[x>pars["high"]]=pr2[2]
          round(pr,round)
        }
      } else fail3=TRUE
    }
    
    
    ### Flat after breakpoint
    if(typ==2 | fail3){
      fct2<-function(x,a,b,brk) ifelse(x <= brk, a+b * x, a+b*brk)
      m <- try(nls(dppm ~ fct2(calmz,a,b,brk),
                   start = list(a = calpars["a"], b = calpars["b"],brk=calpars["brk"]),
                   data=matdf,subset = out),TRUE)
      if(!"try-error"%in%class(m)){
        calpars<-c(coef(m),low=calInt[1])
        calfct<-function(x,pars,round=3){
          pr=ifelse(x <= pars["brk"], pars["a"]+pars["b"]*x, pars["a"]+pars["b"]*pars["brk"])
          pr[x<pars["low"]]=pars["a"]+pars["b"]*pars["low"]
          round(pr,round)
        }
      }
    }
    matdf$pred=calfct(matdf$calmz,calpars)
    matdf$res=matdf$dppm-matdf$pred
    
    return(list(matdf,calpars,calfct))
    
  }
  
  if(is.character(lmzcalib) & length(lmzcalib)==1){
    data(BGIons)
    if(!lmzcalib%in%names(BGIons)) stop('Method not recognised')
    cat("Using predefined background ions for method ",lmzcalib," -> n=",sep="")
    lmzcalib=BGIons[[lmzcalib]]
    cat(length(lmzcalib),"\n")
  }
  
  # dppm=20;bw=2.5;minpts=31;bkpoint=200;calInt=c(80,350)
  
  if("xcmsRaw"%in%class(filein)){
    xRaw=filein
    cat("xRaw file from ",unclass(xRaw@filepath),"\n",sep="")
  }
  if("character"%in%class(filein)){
    cat("Opening ",filein,"\n")
    xRaw=xcmsRaw(filein,includeMSn = FALSE)
  }
  allcal=list()
  for(ical in 1:length(lmzcalib)){
    imz=.GRrawMat(xRaw,mzrange=lmzcalib[ical]*(1+c(-1,1)*dppm*10^-6))
    if(nrow(imz)<minpts) next
    allcal[[ical]]=cbind(imz,cal=ical,cl=.GRsplitOne(imz[,"mz"],bw=bw,ismass = T)$cl)
    # plot(imz[,"rt"],10^6*(1-imz[,"mz"]/ical),col=imz[,"cols"],main=ical,ylim=c(-dppm,dppm))
    # abline(h=0)
  }
  matcal=do.call("rbind",allcal)
  cal2cl=paste(matcal[,"cal"],matcal[,"cl"],sep=";")
  cat(" --found ",length(unique(matcal[,"cal"])),"/",length(lmzcalib)," and ",length(unique(cal2cl))," mass groups\n",sep="")
  matdf=data.frame(do.call("rbind",tapply(1:nrow(matcal),cal2cl,function(x) 
    c(mz=matrixStats::weightedMedian(matcal[x,"mz"],matcal[x,"y"]),n=length(x),matcal[x[1],c("cal" , "cl") ]))))
  matdf$calmz=lmzcalib[matdf$cal]
  matdf$dppm=10^6*(matdf$mz/matdf$calmz-1)
  matdf=matdf[matdf$n>=minpts,]
  matdf=matdf[order(matdf$calmz,matdf$mz),]
  if(retMat) return(matdf)
  
  ### estimate outliers around breakpoint
  matdf$res=0
  l1=which(matdf$calmz<=bkpoint)
  matdf$res[l1]=resid(MASS:::rlm(dppm~calmz,matdf[l1,]))
  l2=which(matdf$calmz>bkpoint)
  matdf$res[l2]=resid(MASS:::rlm(dppm~1,matdf[l2,]))
  matdf$out=(1:nrow(matdf)%in%tapply(1:nrow(matdf),matdf$cal,function(x){
    minres=min(abs(matdf$res[x]))
    x[abs(matdf$res[x])<(2*minres)]}))

  ### fit segments+estimate breakpoints
  mod=segmented(lm(dppm~calmz,data=matdf,subset = out),seg.Z = ~calmz)
  out=car:::outlierTest(lm(resid(mod)~1), cutoff=0.05,n.max=10)
  while(any(out$signif)){
    matdf[names(out$bonf.p),]$out=FALSE
    mod=segmented(lm(dppm~calmz,data=matdf,subset = out),seg.Z = ~calmz)
    out=car:::outlierTest(lm(resid(mod)~1), cutoff=0.05,n.max=10)
  }
  
  ### Regressing around breakpoint
  calpars<-c(a=unname(coef(mod)[1]),b=unname(coef(mod)[2]),c=unname(coef(mod)[3]),brk=unname(mod$psi[2]),low=calInt[1],high=calInt[2])
  calfct<-function(x,pars,round=3){
    pr=pars["a"]+pars["b"]*x
    pr[x>pars["brk"]]=(pr+pars["c"]*(x-pars["brk"]))[x>pars["brk"]]
    pr[x<pars["low"]]=pars["a"]+pars["b"]*pars["low"]
    pr[x>pars["high"]]=pars["a"]+pars["b"]*pars["high"]+pars["c"]*(pars["high"]-pars["brk"])
    round(pr,round)
  }

  matdf$pred=calfct(matdf$calmz,calpars)
  matdf$res=matdf$dppm-matdf$pred
  oldout=rep(FALSE,nrow(matdf))
  
  ## loop over since pb with duplicated mass flagged as outlying
  while(typ>1 & sum(oldout!=matdf$out)){
    oldout=matdf$out
    re=.infct(matdf,calpars,calfct,typ=typ)
    matdf=re[[1]]
    calpars=re[[2]]
    calfct=re[[3]]
    l2chk=unique(intersect(matdf$calmz[matdf$out],matdf$calmz[!matdf$out]))
    if(length(l2chk)>0){
      for(i in l2chk){
      l=which(matdf$calmz==i)
      matdf$out[l[which.min(abs(matdf$res[l]))]]=TRUE
      matdf$out[l[-which.min(abs(matdf$res[l]))]]=FALSE
      }
    }
  }
  
  
  if(doPlot){
    par(mfrow=c(1,2),mar=c(5,4.5,1,.2),cex.main=1)
    xl=pretty(range(c(matdf$calmz,calInt)))
    pred=calfct(min(xl):max(xl),calpars)
    yl=pretty(range(c(matdf$dppm,0,pred),na.rm=T))
    plot(matdf$calmz,matdf$dppm,col=matdf$out+1,ylim=range(yl),axes=F,xlim=range(xl),xlab="m/z",ylab="delta ppm",main=basename(filein))
    abline(h=-1:1,lty=c(2,1,2))
    abline(v=calpars["brk"],col="grey",lwd=2)
    abline(v=calpars[names(calpars)%in%c("low","high")],col="grey",lwd=1,lty=2)
    axis(2,at=yl,las=2,pos=min(xl));axis(1,at=xl,pos=min(yl))
    lines(min(xl):max(xl),pred,col=2,lwd=2)
    
    pred=predict(mgcv:::gam(res~s(calmz),data=matdf,subset=out),newdata = data.frame(calmz=min(xl):max(xl)))
    yl=pretty(range(c(matdf$res[matdf$out],-1:1,pred),na.rm=T))
    plot(matdf$calmz,matdf$res,col=matdf$out+1,ylim=range(yl),axes=F,xlim=range(xl),xlab="m/z",ylab="Residuals",main=basename(filein))
    lines(min(xl):max(xl),pred,col=4,lwd=2)
    abline(h=-1:1,lty=c(2,1,2))
    abline(v=calpars["brk"],col="grey",lwd=2)
    abline(v=calpars[names(calpars)%in%c("low","high")],col="grey",lwd=1,lty=2)
    axis(2,at=yl,las=2,pos=min(xl));axis(1,at=xl,pos=min(yl))
    par(mfrow=c(1,1))
  }
  list(mat=matdf,calfct=calfct,calpars=calpars)
  # mod=scam(dppm~s(calmz,bs="mpi",m=2),family=gaussian(link="identity"),data=matdf[matdf$out,],not.exp=FALSE)
  # modd=scam(dppm~s(calmz,bs="mpd",m=2),family=gaussian(link="identity"),data=matdf[matdf$out,],not.exp=FALSE)
  # if(modd$gcv.ubre<mod$gcv.ubre) mod=modd
  
  # pred=predict(mod,predType = "full-model",newdata = data.frame(calmz=1:ceiling(max(xRaw@mzrange))))
  # names(pred)=1:ceiling(max(xRaw@mzrange))
  # n1=max(calInt[1],ceiling(min(matdf$calmz[matdf$out])))
  # n2=min(calInt[2],floor(max(matdf$calmz[matdf$out])))
  # pred=pred[n1:n2]
#  list(mat=matdf,pred=pred,coef=c(coef(mod),n1=n1,n1dppm=unname(pred[1]),n2=n2,n2dppm=unname(rev(pred)[1])))
}

######################################################################################
## Single point MZ calib over scan/mz
#' get the xcmsRaw lighter
#'
#' @name .MGcorrRawPPM
#' @rdname .MGcorrRawPPM-methods
#' @import xcms
#' @export
setGeneric(name=".MGcorrRawPPM",def=function(object,calfct,calpars) standardGeneric(".MGcorrRawPPM"))


#' @rdname .MGcorrRawPPM-methods
#' @aliases .MGcorrRawPPM
setMethod(".MGcorrRawPPM","xcmsRaw",
          .MGcorrRawPPM<-function(object,calfct,calpars){
            ## based on calibration function/parameeters

           if("character"%in%class(object)){
              cat("Opening ",object,"\n")
              object=xcmsRaw(object,includeMSn = FALSE)
            }
            
            scsten=cbind(object@scanindex+1,c(object@scanindex[-1],length(object@env$intensity)))
            rownames(scsten)=1:length(object@scanindex)

            newmz=object@env$mz
            for(ix in rownames(scsten)){
              oldy<-xcms:::getScan(object,as.numeric(ix))[,1]
              ppmy<- 1-calfct(oldy,calpars) *10^-6
              #              ppmy<- 1-addppm*10^-6
              newy=round(oldy*ppmy,6)
              newmz[scsten[ix,1]:scsten[ix,2]]=newy
            }
            print(summary((1-newmz/object@env$mz)*10^6))
            
            ob<-new("xcmsRaw")
            ob@env <- new.env(parent = .GlobalEnv)
            ob@env$mz<-as.numeric(newmz)
            ob@env$intensity<-object@env$intensity
            ob@scanindex<-object@scanindex
            ob@scantime<-object@scantime
            
            ob@acquisitionNum<-1:length(ob@scanindex)
            ob@filepath<-object@filepath
            ob@mzrange<-range(ob@env$mz)
            ob@profmethod<-object@profmethod
            ob@tic<-object@tic
            ob@profparam<-list()
            ob<-xcms:::remakeTIC(ob)
            return(ob)
          }
          
)
######################################################################################
