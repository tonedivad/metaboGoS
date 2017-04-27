######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
### Serious 3.4.1 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' Interface to Sirius
#'
#' Interface to Sirus
#'
#' @param precmz Precursor M/Z
#' @param ms1spec MS1 spectra
#' @param ms2spec MS/MS spectra
#' @param ionProd Ionisation product
#' @param ppmPrecMF ppmPrecMF
#' @param maxMF Maximum number of MF to be returned
#' @param limint Minimum intensity perc.
#' @param listmf Annotation of a single list of mol formula
#' @param retAnn Add annotation in the original MS/MS
#' @param exec Sirius executable
#' @param tempdir Directory containing the output
#' @param basenam result file name
#' @param clean  clean-up the mess
#' @return List of stuff
#' @export
getMF.sirius<-function(precmz=NULL,ms1spec=NULL,ms2spec=NULL,ionProd="[M+H]1+",ppmPrecMF=11,maxMF=NA,limint=10^-5,
                      listmf=c(),retAnn=FALSE,
                     exec="/media/D01/Metabo/MetSoft/sirius3.4.1/bin/sirius3",
                     tempdir="./",basenam=paste0("tmp",round(runif(1)*1000)),clean=TRUE){
  
  ### Set tempfiles
  ifile=paste0(tempdir,"/",basenam,".ms")
  
  allre=allann=list()
  rematch= data.frame("PrecMZ"=double(),"Mass"=double(),"PPM"=double(),
                      "MF"=character(),"Adduct"=character(),"Score"=double(),"
                      ScoreMS2"=double(),"ScoreIso"=double(),"NPeakMS2"=integer(),
                      "PMatch"=double(),"NPeakIso"=integer())
  ann=NULL
  
      cmpdnam=sprintf('MyPrec%d',round(precmz*10000))
      cat(">compound ", cmpdnam,"\n",sep="",file = ifile,append = F)
      cat(">parentmass ", precmz,"\n",sep="",file = ifile,append = T)
      cat(">ionization ", gsub("(.*)[0-9]([\\+-])$","\\1\\2",ionProd),"\n",sep="",file = ifile,append = T)
      
      ## Add MS1
      if(!is.null(ms1spec)){
        cat("\n>ms1\n",sep="",file = ifile,append = T)
        for(i in which((ms1spec[,"y"]/max(ms1spec[,"y"]))>limint)) cat(round(ms1spec[i,"mz"],5)," ",round(ms1spec[i,"y"],5),"\n",sep="",file = ifile,append = T)
      }
      
      ## Add MS2
      if(!is.null(ms2spec)){
        cat("\n>ms2\n",sep="",file = ifile,append = T)
      for(i in which((ms2spec[,"y"]/max(ms2spec[,"y"]))>limint)) cat(round(ms2spec[i,"mz"],5)," ",round(ms2spec[i,"y"],5),"\n",sep="",file = ifile,append = T)
      }
      
      args0=c(paste0("-p ","orbitrap"),"-s score",paste0("--ppm-max ",ppmPrecMF))
      if(length(listmf)==0)   args0=c(args0,paste0("-c ",ifelse(is.na(maxMF),ceiling(precmz),maxMF)))
      if(length(listmf)>0 | retAnn){
        retAnn=TRUE
        outdir=paste0(tempdir,"/",basenam,"anndir")
        if(!dir.exists(outdir)) dir.create(outdir)
        if(length(listmf)>0) args0=c(paste0("-f ",paste(listmf,collapse = " ")),args0)
        args0=c(args0,paste0("-o ",outdir))
      }
      args=c(args0,ifile)
     print(paste(args,collapse = " "))
      ## parse results
      re<-system2(exec,args=args,stdout = T,stderr = F)
      re=re[grep("score",re)]
      if(length(re)==0){
        if(clean) file.remove(ifile)
        if(clean & retAnn) unlink(outdir,recursive = T)
        return(list(Match=rematch,Ann=ann))
      }
      rematch=data.frame(do.call("rbind",lapply(strsplit(re,"\t"),function(x) gsub(".*: ","",x))))
      names(rematch)=c("MF","Adduct","Score","ScoreMS2","ScoreIso","NPeakMS2","PMatch","NPeakIso")
      for(i in 3:ncol(rematch)) rematch[,i]=as.numeric(gsub(" .*","",rematch[,i]))
      rematch$MF=gsub(".* ","",rematch$MF)
      rematch$Adduct=ionProd
      newmass=round(sapply(rematch$MF,function(x) rcdk:::get.formula(x)@mass),6)
      rematch=data.frame("PrecMZ"=precmz,"Mass"=newmass,"PPM"=round((newmass/precmz-1)*10^6,4),rematch)
    
       if(retAnn){
        
        lfiles=list.files(paste0(outdir,"/1_",basenam,"_",cmpdnam,"/spectra"),full.names = T)
        df=lapply(lfiles,function(ifi){
          df=data.frame(read.table(ifi,sep="\t",header=T))
          df=df[match(1:nrow(ms2spec),apply(abs(outer(ms2spec[,"mz"],df$mz,"-")),2,which.min)),c("explanation","exactmass")]
          df})
        names(df)=gsub(".*_([A-Z0-9]+)\\.ms$","\\1",lfiles)
        df=df[match(rematch$MF,names(df))]
        lnames=names(df)
        df=do.call("cbind",df)
        colnames(df)=paste0(rep(lnames,each=2),c(".MF",".mz"))
        ann=data.frame(ms2spec,df)
        ann=ann[rowSums(is.na(df))!=ncol(df),]
        rownames(ann)=NULL
      }
      if(clean) file.remove(ifile)
      if(clean & retAnn) unlink(outdir,recursive = T)
      
  invisible(list(Match=rematch,Ann=ann))
}



######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
### 7GR 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' Interface to the 7GR executable
#'
#' Interface to 7GR executable
#'
#' @param precmz Precursor m/z
#' @param lIP Vector of ionisation products
#' @param precppm Accuracy in ppm
#' @param precmtol Tolerence in mDa
#' @param ipdb IP database
#' @param exec 7GR executable
#' @param chkdb database of molform
#' @import rcdk
#' @return List of stuff
#' @export
getMF.7GR<-function(precmz,lIP="[M+H]1+",precppm=11,precmtol=10,ipdb=NULL,exec='/media/D01/Metabo/MetSoft/HR/myhr',chkdb=NULL){
  
  # lIP="[M+H]1+";precppm=11;precmtol=10;ipdb=IPDB;exec='/media/D01/Metabo/MetSoft/HR/myhr';chkdb="MolFormDF"
  .inHR<-function(imz,mtolprec,exec='myhr'){
    
    args=c(paste0("-m ",imz),paste0("-t ",mtolprec))
    system.time(re<-system2(exec,args=args,stdout = T))
    re=strsplit(re,"\t")
    re=re[sapply(re,length)==7]
    if(length(re)==0) return(NULL)
    re=do.call("rbind",re)[,2:7,drop=F]
    re=data.frame(do.call("rbind",lapply(re[,1],function(x){x=get.formula(x);data.frame(Query=imz,MF=x@string,Mass=x@mass)})),RDB=re[,2])
    re$mTol=(re$Mass-imz)*1000
    re$Dppm=10^3*(re$mTol)/imz
    re=re[order(abs(re$mTol)),]
    invisible(re)
  }
  
  if (is.null(ipdb)) {
    data(IPDB)
    ipdb = IPDB
  }
  
  
  ladd = unique(which(ipdb$Id %in% lIP | ipdb$Name %in% lIP | ipdb$Set %in% lIP))
  
  if (length(ladd) == 0) 
    stop(paste("IPs invalid - must be in :\n", "  *positive mode: ", 
               paste(IPDB$Name[IPDB$Charge > 0], collapse = " "), 
               "\n", "  *negative mode: ", paste(IPDB$Name[IPDB$Charge < 
                                                             0], collapse = " "), "\n", sep = ""))
  
  
  m2match=do.call("rbind",lapply(ladd, function(i) cbind(i,precmz,(abs(ipdb[i, ]$Charge) * precmz - ipdb[i, ]$adMass)/ipdb[i, ]$xM)))
  if(!is.null(precppm)) m2match=cbind(m2match,precmtol=m2match[,3]*precppm*10^-3) else m2match=cbind(m2match,precmtol=precppm)
  
  allre=list()
  for(i in 1:nrow(m2match)){
    re=.inHR(imz = m2match[i,3],mtolprec = m2match[i,4],exec = exec)
    if(is.null(re)) next
    re$Adduct=ipdb[m2match[i,1],]$Name
    re$precMZ=m2match[i,2]
    allre[[i]]=re[,c(  "Adduct", "precMZ","MF","Mass","RDB","mTol","Dppm")]
  }
  res=do.call("rbind",allre)
  rownames(res)=NULL
  if(!is.null(chkdb)){
    eval(parse(text=paste0("data(",chkdb,")")))
  res$InDB=res$MF%in%MolFormDF$MF
  if(any(res$Adduct=="[M]1+")) res$InDB[res$Adduct=="[M]1+"]=res$MF[res$Adduct=="[M]1+"]%in%MolFormDF$MF[MolFormDF$Charge==1]
  if(any(res$Adduct=="[M]1-")) res$InDB[res$Adduct=="[M]1-"]=res$MF[res$Adduct=="[M]1-"]%in%MolFormDF$MF[MolFormDF$Charge== -1]
  res=res[order(!res$InDB),,drop=F]
  }
  return(res)
}


