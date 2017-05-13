### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name makeMetaboSet
#' @title makeMetaboSet
#' 
#' stuff
#'
#' @param matpks peak matrix
#' @param metainfos metainfos
#' @return metaboSet
#' @export
makeMetaboSet<-function(matpks,metainfos){
  imeth=unique(metainfos$Method)
  lusampid=metainfos$Sid
  # metainfos=metaData[which(metaData$Sid%in%lusampid),c("Sid","sType","InjOrder")]
  # fileinfos=metaData[which(metaData$Sid%in%lusampid),c("Sid","dirName","fileName","completionTime")]
  # lusampid=metainfos$Sid
  
  ### define new metabolites
  newpkid=paste(matpks$RoiId,matpks$Pk,sep=";;")
  l=which(matpks$InDeco==1)
  newpkids=do.call("rbind",tapply(l,newpkid[l],function(x) 
    data.frame(PkId=newpkid[x][1],mz=matrixStats::weightedMedian(matpks[x,"mzmed"],matpks[x,"int.sm"]),
               rt=matrixStats::weightedMedian(matpks[x,"rt"],matpks[x,"int.sm"]))))
  newpkids$mz=round(newpkids$mz,6)
  newpkids$rt=round(newpkids$rt,3)
  
  ### Prep data
  pkid2sid=tapply(1:length(newpkid),newpkid,function(x) x[match(lusampid,matpks$Sid[x])])
  newpkids=newpkids[match(names(pkid2sid),newpkids$PkId),]
  newpkids$metnam=sprintf("%.5f@%.3f-%s",newpkids$mz,newpkids$rt,imeth)
  ldups=names(which(table(newpkids$metnam)>1))
 if(length(ldups)>0) for(idup in ldups){
  l=which(newpkids$metnam==idup)
  newpkids$metnam[l]=sprintf("%.5f-%d@%.3f-%s",newpkids$mz[l],1:length(l),newpkids$rt[l],imeth)
}
  lvars=names(matpks)[!names(matpks)%in%c("RoiId","Sid","Pk","InDeco")]
  
  alldata=list()
  for(ivar in lvars){
    m=do.call("cbind",lapply(pkid2sid,function(x) round(matpks[x,ivar],6)))
    dimnames(m)=list(lusampid,newpkids$metnam)
    alldata[[ivar]]=m
  }
  if("rt"%in%names(alldata) & !"RT"%in%names(alldata)) names(alldata)[names(alldata)=="rt"]="RT"
  if("mz"%in%names(alldata) & !"MZ"%in%names(alldata)) names(alldata)[names(alldata)=="mz"]="MZ"
  if("int.ap"%in%names(alldata) & !"Height"%in%names(alldata)) names(alldata)[names(alldata)=="int.ap"]="Height"
  if("area"%in%names(alldata) & !"Area"%in%names(alldata)) names(alldata)[names(alldata)=="area"]="Area"
  ### Make Annot
  annot=data.frame(Analyte=newpkids$metnam,MetName=NA,LevelAnnot=4,Method=imeth,IsSTD=FALSE,IsISO=FALSE,RT=newpkids$rt,MZ=newpkids$mz,PkId=newpkids$PkId)
  rownames(annot)=annot$Analyte
  ### Make meta
  metadf=metainfos[,!names(metainfos)%in%c(  "dirName","fileName","completionTime" ,"Method")]
  
  filedf=metainfos[,names(metainfos)%in%c(  "Sid","sType")]
  filedf$File=paste(metainfos$dirName,metainfos$fileName,sep="/")
  filedf$Date=metainfos$completionTime
  if(!is.null(metainfos$Batch)) filedf$Batch=metainfos$Batch
  
  filedf$InjOrder=metadf$InjOrder=order(order(filedf$Date))
  
  allmat=list(Method=imeth,Sid=lusampid,Analyte=annot$Analyte,Annot=annot,Meta=metadf,File=filedf,Data=alldata)
  class(allmat)=append(class(allmat),"metaboSet")
  invisible(allmat)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' @name mmapMS2toMetaboSet
#' @title mmapMS2toMetaboSet
#' 
#' stuff
#'
#' @param cmb metaboSet object
#' @param fileres list DDA files
#' @param dmz delta mz 
#' @param dppm delta ppm
#' @param rtwin rtwin 
#' @param outfile file to be exported
#' @return metaboSet
#' @export

mapMS2toMetaboSet<-function(cmb,fileres,dmz=0.001,dppm=2,rtwin=3/60,outfile=NULL){
  
  # dmz=0.001;dppm=2;rtwin=3/60
  
  pk2m=data.frame(Ana=cmb$Analyte,
                  rtmin=apply(cmb$Data$rtmin[names(fileres),],2,min,na.rm=T),
                  rtmax=apply(cmb$Data$rtmax[names(fileres),],2,max,na.rm=T),mzmed=cmb$Annot$MZ,
                  rtmin0=NA,rtmax0=NA)
  rownames(pk2m)=pk2m$Ana
  
  cmbPrec2Ana=cmbMS2Data=list()
  for(isid in names(fileres)){
    cat("Processing ",isid,sep="")
    load(fileres[isid])
    pk2m$rtmin0=cmb$Data$rtmin[isid,]
    pk2m$rtmax0=cmb$Data$rtmax[isid,]
    
    ll2sp=vector("list",nrow(pk2m))
    names(ll2sp)=rownames(pk2m)
    
    ms2p=obj$MS2toPrec
    for(i in rownames(pk2m)){
      ddrt=(ms2p$rt-pk2m[i,"rtmin"])>=(-rtwin) & (ms2p$rt-pk2m[i,"rtmax"])<=(rtwin)
      ddmz=(abs(pk2m[i,]$mzmed-ms2p$mz)<dmz | (abs(1-pk2m[i,]$mzmed/ms2p$mz)<=dppm*10^-6))
      l=which(ddmz & ddrt)
      if(length(l)==0) next
      inpk=(ms2p$rt[l]>=pk2m[i,]$rtmin0 & ms2p$rt[l]<=pk2m[i,]$rtmax0)
      inpkg=(ms2p$rt[l]>=pk2m[i,]$rtmin & ms2p$rt[l]<=pk2m[i,]$rtmax)
      case=any(inpk)*1+any(inpkg)*1
      case=ifelse(is.na(case),0,case)
      if(case==1)  l=l[which(inpkg)]
      if(case==2)  l=l[which(inpk)]
      ll2sp[[i]]=data.frame(Analyte=i,Case=case,i=tapply(l,ms2p$SpId[l],function(x) x[which.max(ms2p$pint[x])]))
    }
    nMS2pk=do.call("rbind",ll2sp)
    nMS2pk=cbind(nMS2pk,ms2p[nMS2pk$i,])
    cmbPrec2Ana[[isid]]=nMS2pk
    cat(": ",sum(unique(obj$Sp2Pk$SpId)%in%nMS2pk$SpId)," out of ",length(unique(unique(obj$Sp2Pk$SpId))),"\n",sep="")
    cmbMS2Data=c(cmbMS2Data,obj$MS2Data[which(names(obj$MS2Data)%in%nMS2pk$SpId)])
  }
  cmbPrec2Ana=do.call("rbind",cmbPrec2Ana)
  rownames(cmbPrec2Ana)=NULL
  cat(" --> ",length(unique(cmbPrec2Ana$Analyte))," peaks, ", nrow(cmbPrec2Ana)," prec, ",length(unique(cmbPrec2Ana$SpId))," MS/MS\n",sep="")
  
  if(!is.null(outfile)) save(file=outfile,cmb,cmbPrec2Ana,cmbMS2Data)
  
}