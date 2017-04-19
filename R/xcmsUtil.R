#' get the xcmsRaw lighter
#'
#' @name eicXcmsRaw
#' @rdname eicXcmsRaw-methods
#' @import xcms
#' @export
setGeneric(name="eicXcmsRaw",def=function(object,eicmat,inthr) standardGeneric("eicXcmsRaw"))


#' @rdname eicXcmsRaw-methods
#' @aliases eicXcmsRaw
setMethod("eicXcmsRaw","xcmsRaw",
          eicXcmsRaw<-function(object,eicmat,inthr=0){
            
            
            ascans=do.call("rbind",lapply(1:nrow(eicmat),function(i)
              .GRrawMat(object,mzrange = eicmat[i,c("mzmin","mzmax")],rtrange = eicmat[i,c("rtmin","rtmax")]*60)[,1:4,drop=F]))
            ascans=ascans[order(ascans[,1],ascans[,2]),,drop=FALSE]
            ascans=ascans[!duplicated(ascans[,4]),,drop=FALSE]
            ascans=ascans[which(ascans[,"y"]>inthr),,drop=FALSE]
            scrange=min(ascans[,1]):max(ascans[,1])
            sc2time=object@scantime[scrange]
            
            tabsc=table(factor(ascans[,1],levels = scrange))
            scanidx=c(0,cumsum(tabsc))
            scanidx=as.integer(scanidx[-length(scanidx)])
            
            newtic=round(as.numeric(tapply(ascans[,"y"],ascans[,"scan"],sum)[as.character(scrange)]),1)
            newtic[is.na(newtic)]=0
            
            
            ob<-new("xcmsRaw")
            ob@env <- new.env(parent = .GlobalEnv)
            ob@env$mz<-as.numeric(round(ascans[,"mz"],6))
            ob@env$intensity<-as.numeric(round(ascans[,"y"],1))
            ob@scanindex<-scanidx
            ob@scantime<-sc2time
            
            ob@acquisitionNum<-1:length(ob@scanindex)
            ob@filepath<-object@filepath
            ob@mzrange<-range(ob@env$mz)
            ob@profmethod<-object@profmethod
            ob@tic<-newtic
            ob@profparam<-list()
            ob<-xcms:::remakeTIC(ob)
            
            return(ob)
          }
          
)

######################################################################################
#' @title Remove artefact in xcmsRaw object
#' get the xcmsRaw lighter
#'
#' @name satremXcmsRaw
#' @rdname satremXcmsRaw-methods
#' @param object xcmsRaw object
#' @param reso resolution
#' @param typ 1/2 2 correrspond to the extended model
#' @param inthr intensity threshold
#' @param rtlim min/max retention time in min.
#' @import xcms
#' @export
setGeneric(name="satremXcmsRaw",def=function(object,reso,typ,inthr,rtlim) standardGeneric("satremXcmsRaw"))


#' @rdname satremXcmsRaw-methods
#' @aliases satremXcmsRaw
setMethod("satremXcmsRaw","xcmsRaw",
          satremXcmsRaw<-function(object,reso=70000,typ=2,inthr=0,rtlim=c(-Inf,Inf)){
            
            scsten=cbind(object@scanindex+1,c(object@scanindex[-1],length(object@env$intensity)))
            rownames(scsten)=1:length(object@scanindex)
            
            newdat=list()
            nrm1=nrm2=0
            lsc=which((object@scantime/60)>=rtlim[1] & (object@scantime/60)<=rtlim[2])
            cat(" -- processing",length(lsc),"scans\n")
            for(ix in lsc){
              re<-xcms:::getScan(object,as.numeric(ix))
              l2excl=which(re[,2]<inthr)
              nrm1=nrm1+length(l2excl)
              if(length(l2excl)>0) re=re[-l2excl,,drop=F]
              l2excl=satremoval(re,reso,typ)
              nrm2=nrm2+length(l2excl)
              if(length(l2excl)>0) re=re[-l2excl,,drop=F]
              if(nrow(re)==0) next
              newdat[[ix]]=cbind(re,scan=ix)
            }
            ntot=diff(range(scsten[as.character(lsc),]))+1
            cat(" -- removed",nrm1,"/",nrm2,"out of",ntot,"duples - avg.=",round(100*nrm1/ntot,2),"/",round(100*nrm2/ntot,2),"perc. per scan")
            newdat=do.call("rbind",newdat)
            newdat=newdat[order(newdat[,"scan"],newdat[,1]),,drop=F]
           scrange=min(newdat[,"scan"]):max(newdat[,"scan"])
            sc2time=object@scantime[scrange]
            
            tabsc=table(factor(newdat[,"scan"],levels = scrange))
            scanidx=c(0,cumsum(tabsc))
            scanidx=as.integer(scanidx[-length(scanidx)])
            
            newtic=round(as.numeric(tapply(newdat[,"intensity"],newdat[,"scan"],sum)[as.character(scrange)]),1)
            newtic[is.na(newtic)]=0
            
            
            ob<-new("xcmsRaw")
            ob@env <- new.env(parent = .GlobalEnv)
            ob@env$mz<-as.numeric(round(newdat[,"mz"],6))
            ob@env$intensity<-as.numeric(round(newdat[,"intensity"],1))
            ob@scanindex<-scanidx
            ob@scantime<-sc2time
            
            ob@acquisitionNum<-1:length(ob@scanindex)
            ob@filepath<-object@filepath
            ob@mzrange<-range(ob@env$mz)
            ob@profmethod<-object@profmethod
            ob@tic<-newtic
            ob@profparam<-list()
            ob<-xcms:::remakeTIC(ob)
            
            return(ob)
          }
          
)
######################################################################################
#' Artfact removal in duples
#' @param re matrix of mz/int
#' @param reso mass resolution
#' @param typ 1/2 extended
#' @return ids of entries to be remove
#' @export
satremoval<-function(re,reso,typ=2){
  #    print(reso)
  if(is.na(reso)) return(c())
  m=re[,1]
  v=re[,2]
  
  if(typ==1) unique(unlist(lapply(1:length(m),function(i){
    m0=m[i]
    gamma=m0/reso
    v1=v[i]*(gamma^2)/(gamma^2+(m-m0)^2)
    which(v<(0.99*v1))})))
  
  if(typ==2) unique(unlist(lapply(1:length(m),function(i){
    m0=m[i]
    gamma=m0/reso
    v1=v[i]*(gamma^2)/(gamma^2+(m-m0)^2)
    gamma2=m0/(reso*0.05)
    v2=0.05*v[i]*(gamma2^2)/(gamma2^2+(m-m0)^2)
    # l=which(v1<(0.05*v[i]))
    # v1[l]=v2[l]
    which(v<(0.99*v1) | v<(0.99*v2))})))
  
  
}


