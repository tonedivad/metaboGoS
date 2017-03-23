#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#' set-up Tremolo input
#'
#' Parse Tremolo results
#'
#' @param what Object containg spectra
#' @param lsp2exp List of spectra
#' @param file2exp where to write the MGF formatted spectra
#' @param perclim Minimum intensity in percentage
#' @param nmin Minimum number of peaks higher than perclim
#' @param exec Conversion executable
#' @return List of stuff
#' @export
getTremololist<-function(what,lsp2exp,file2exp=NULL,perclim=0.1,nmin=5,exec='/media/D01/Metabo/MetSoft/Tremolo_8_21_2013/convert'){
  
  if(!"ddaSet"%in%class(what)) stop("Error: not a ddaSet object")
  
  toexp=alladd=list()
  iix=0
  for(i in 1:length(lsp2exp)){
    isp=lsp2exp[i]
    pol=what$MS2Infos[isp,]$Polarity
    m=what$MS2Data[[isp]][,c("mz","y")]
    m[,"y"]=round(100*m[,"y"]/max(m[,"y"]),4)
    m=m[m[,2]>=perclim,,drop=F]
    if(nrow(m)<nmin) next
    toadd=what$MS2Infos[which(what$MS2ClInfos$SpId==isp),c("PrecMZ","RT","PrecInt","SpId")]
    if(!is.null(what$MS2ClInfos)) toadd=what$MS2ClInfos[which(what$MS2ClInfos$SpId==isp),c("MZ","RT","He","SpId")]
    for(j in 1:nrow(toadd)){
      iix=iix+1
      alladd[[iix]]=cbind(SCANS=iix,toadd[j,])
      toexp[[iix]]=c('BEGIN IONS',sprintf("PEPMASS=%.5f",toadd[j,1] ),'SEQ=*..*',
                     sprintf("CHARGE=%s",ifelse(pol=="Negative","1-","1+")),
                     sprintf("IONMODE=%s", pol),
                     sprintf("TITLE= MS/MS at %.2f int=%e",  toadd[j,2],round(toadd[j,3])),
                     sprintf("SCANS=%d",iix),
                     apply(m,1,function(x) sprintf("%.5f %.3f",x[1],x[2])),
                     'END IONS','')
    }
  }
  alladd=do.call("rbind",alladd)
  rownames(alladd)=NULL
  if(!is.null(file2exp) & length(toexp)>0){
    write.table(alladd,file=paste0(file2exp,".pkinfo"),sep="\t",row.names = F,col.names = F)
    cat(unlist(toexp),file=paste0(file2exp,".mgf"),sep="\n")
    system(paste0(exec," ",file2exp,".mgf"))
  }
  invisible(list(alladd,toexp))
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#' read Tremolo output
#'
#' Parse Tremolo results
#'
#' @param root Directory to store the results
#' @param speclib Spectra lib in .mgf format
#' @param prectol Precursor tolerance
#' @param score Minimum score
#' @param topk Only keep the top matches
#' @param exec Tremolo execuytable
#' @return Data frame
#' @export
runTremolo<-function(root,speclib,prectol=0.1,score=0.4,topk=10,exec='/media/D01/Metabo/MetSoft/Tremolo_8_21_2013/main_execmodule'){
  
  
  # root="./tmp"
  # speclib="../ConvSMi/AnnotSpecBA.mgf"
  # prectol=.2
  # score=0.4
  # topk=50
  
  cat(paste0("EXISTING_LIBRARY_MGF=",normalizePath(speclib),"\n\n","searchspectra=",root,".pklbin\n\n","RESULTS_DIR=",root,".out0\n\n",
             'tolerance.PM_tolerance=',prectol,"\n\nsearch_decoy=0\n\nSCORE_THRESHOLD=",score,"\nTOP_K_RESULTS=",topk,"\n\n",
             'NODEIDX=0\nNODECOUNT=1\n\nSEARCHABUNDANCE=0\nSLGFLOADMODE=1'),file=paste0(root,".params"))
  
  system(paste0(exec,' ExecSpectralLibrarySearch ',paste0(root,".params")))
  system(paste0('cut -f 1,6,11,16,21,25 ',root,".out0 > ",root,".out"))
  
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#' read Tremolo output
#'
#' Parse Tremolo results
#'
#' @param root Directory to clean
#' @return Data frame
#' @export
readTremolo<-function(root){
  
  outdat=do.call("rbind",strsplit(scan(paste0(root,".out"),what="raw",sep="\n"),"\t"))
  colnames(outdat)=c("Sp",outdat[1,-1])
  outdat=outdat[-1,]
  if(nrow(outdat)==0) return(NULL)
  outdat=data.frame(outdat)
  outdat$MQScore=as.numeric(outdat$MQScore)
  outdat$mzErrorPPM=round(as.numeric(outdat$mzErrorPPM),3)
  outdat$dbIndex=as.numeric(outdat$dbIndex)
  outdat$CmpdEntry=sapply(strsplit(outdat$CompoundName,":"),function(x) x[1])
  outdat$Energy=sapply(strsplit(outdat$CompoundName,":"),function(x) x[2])
  #outdat=outdat[order(outdat$Sp,outdat$LibSearchSharedPeaks,outdat$MQScore),]
  l=which(!duplicated(outdat[,c("Sp","CmpdEntry",'LibSearchSharedPeaks')]))
  outdat=outdat[l,]
  outdat$CmpdEntry2=paste(outdat$CmpdEntry,outdat$Sp,sep=";")
  
  nid=do.call("rbind",tapply(1:nrow(outdat),outdat$CmpdEntry2,function(l) l[c(which.max(outdat$MQScore[l]),which.max(outdat$LibSearchSharedPeaks[l]))]))
  outdat=data.frame(Sp=outdat$Sp[nid[,1]],Cmpd=outdat$CmpdEntry[nid[,1]],DPPM=outdat$mzErrorPPM[nid[,1]],
                    MQScore.1=outdat$MQScore[nid[,1]],
                    NMatch.1=outdat$LibSearchSharedPeaks[nid[,1]],
                    Energy.1=outdat$Energy[nid[,1]],
                    MQScore.2=outdat$MQScore[nid[,2]],
                    NMatch.2=outdat$LibSearchSharedPeaks[nid[,2]],
                    Energy.2=outdat$Energy[nid[,2]])
  outdat$Sp=as.numeric(outdat$Sp)
  outdat=outdat[order(outdat$Sp,outdat$MQScore.1),]
  pkinfos=read.table(paste0(root,".pkinfo"),sep="\t")
  pkinfos=pkinfos[match(outdat[,1],pkinfos[,1]),]
  outdat$Sp=pkinfos[,5]
  outdat$PrecMZ=pkinfos[,2]
  outdat$RT=pkinfos[,3]
  rownames(outdat)=NULL
  return(outdat)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#' Clean up Tremolo mess
#'
#' Clean up Tremolo mess
#'
#' @param root Directory to clean
#' @param removeResfile Should the result file be removed
#' @return Nothing
#' @export
cleanTremolo<-function(root,removeResfile=TRUE){
  if(resFile) l2rm=paste0(root,c( ".mgf",".out",".out0",".params",".pkinfo",".pklbin")) else  l2rm=paste0(root,c( ".mgf",".out0",".params",".pklbin"))
  file.remove(l2rm)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
