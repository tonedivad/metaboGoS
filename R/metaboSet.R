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
  annot=data.frame(Analyte=newpkids$metnam,MetName=NA,LevelAnnot=4,Method=imeth,IsSTD=FALSE,IsISO=FALSE,RT=newpkids$rt,MZ=newpkids$mz)
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

