#' get the xcmsRaw lighter
#'
#' @name slimXcmsRaw
#' @rdname slimXcmsRaw-methods
#' @import xcms
#' @export
setGeneric(name="slimXcmsRaw",def=function(object,eicmat,inthr) standardGeneric("slimXcmsRaw"))


#' @rdname slimXcmsRaw-methods
#' @aliases slimXcmsRaw
setMethod("slimXcmsRaw","xcmsRaw",
          slimXcmsRaw<-function(object,eicmat,inthr=0){
            
            
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



