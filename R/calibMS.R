#' @name calibMS1
#' 
#' @title Convert CFM file to MGF
#' some stufff
#'
#' @param filein Path to the input file
#' @param lmzcalib Calibration masses
#' @param dppm Mass error
#' @param bw Bandwidth to group masses
#' @param minpts Minimum of points
#' @param bkpoint breakpoint mass to calculate robust regression
#' @param calInt calibration interval to compute 
#' @import DoseFinding
#' @return xcms
#' @export
calibMS1<-function(filein,lmzcalib,dppm=20,bw=2.5,minpts=11,bkpoint=median(calInt),calInt=c(80,350)){
  
  # dppm=20;bw=2.5;minpts=11;bkpoint=200;calInt=c(80,350)
  cat("Opening ",filein,"\n")
  xRaw=xcmsRaw(filein,includeMSn = FALSE)
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
  
  ### remove duplicated based on robust regression before/after break point
  matdf$res=0
  l1=which(matdf$calmz<=bkpoint)
  matdf$res[l1]=resid(MASS:::rlm(dppm~calmz,matdf[l1,]))
  l2=which(matdf$calmz>bkpoint)
  matdf$res[l2]=resid(MASS:::rlm(dppm~1,matdf[l2,]))
  matdf$out=(1:nrow(matdf)%in%tapply(1:nrow(matdf),matdf$cal,function(x){
    minres=min(abs(matdf$res[x]))
    x[abs(matdf$res[x])<(2*minres)]}))
  plot(matdf$calmz,matdf$dppm,col=matdf$out+1)
  
  ###  gam  
  mod=DoseFinding:::fitMod(calmz,dppm,matdf[matdf$out,],model="sigEmax",bnds = rbind(calInt,range(matdf$dppm[matdf$out])*1.2))
  #mod=DoseFinding:::fitMod(calmz2,dppm,matdf[matdf$out,],model="emax")
  
  
  
  # mod=scam(dppm~s(calmz,bs="mpi",m=2),family=gaussian(link="identity"),data=matdf[matdf$out,],not.exp=FALSE)
  # modd=scam(dppm~s(calmz,bs="mpd",m=2),family=gaussian(link="identity"),data=matdf[matdf$out,],not.exp=FALSE)
  # if(modd$gcv.ubre<mod$gcv.ubre) mod=modd
  
  pred=predict(mod,predType = "full-model",newdata = data.frame(calmz=1:ceiling(max(xRaw@mzrange))))
  names(pred)=1:ceiling(max(xRaw@mzrange))
  n1=max(calInt[1],ceiling(min(matdf$calmz[matdf$out])))
  n2=min(calInt[2],floor(max(matdf$calmz[matdf$out])))
  pred=pred[n1:n2]
  list(mat=matdf,pred=pred,coef=c(coef(mod),n1=n1,n1dppm=unname(pred[1]),n2=n2,n2dppm=unname(rev(pred)[1])))
}

######################################################################################
## Single point MZ calib over scan/mz
#' get the xcmsRaw lighter
#'
#' @name .MGcorrRawPPM
#' @rdname .MGcorrRawPPM-methods
#' @import xcms
#' @export
setGeneric(name=".MGcorrRawPPM",def=function(object,addppm) standardGeneric(".MGcorrRawPPM"))


#' @rdname .MGcorrRawPPM-methods
#' @aliases .MGcorrRawPPM
setMethod(".MGcorrRawPPM","xcmsRaw",
          .MGcorrRawPPM<-function(object,addppm){
            ## single point calibration
            ## if vector, names must match the correspongi m/z
            
            #  xRaw2=xRaw
            scsten=cbind(object@scanindex+1,c(object@scanindex[-1],length(object@env$intensity)))
            rownames(scsten)=1:length(object@scanindex)
            lmz=floor(object@mzrange[1]):ceiling(object@mzrange[2])
            
            ### Fix add ppm
            if(length(addppm)==1) addppm=rep(round(addppm,3),length(lmz)) else   addppm=round(approx(as.numeric(names(addppm)),addppm,lmz,yleft = addppm[1],yright = rev(addppm)[1])$y,3)
            
            newmz=object@env$mz
            for(ix in rownames(scsten)){
              oldy<-xcms:::getScan(object,as.numeric(ix))[,1]
              ppmy<- 1-approx(lmz,addppm,oldy)$y *10^-6
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
