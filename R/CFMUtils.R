#' Convert CFM file to MGF
#'
#' This function loads a
#'
#' @param ifile Path to the input file
#' @param itmp Path to the input file
#' @param minpk Path to the input file
#' @return List of list
#' @export
Cfm2Mgf<-function(ifile,itmp,minpk=3){
  
  .infct<-function(xtmp,deltaMZ=0.001,whichmz="X1",whichint="X2",keepMax=FALSE){
    # deltaMZ=0.001;whichmz="X1";whichint="X2";keepMax=FALSE;xtmp=tmp
    
    if(nrow(xtmp)<2) return(xtmp)
    xtmp=xtmp[order(xtmp[,whichmz],xtmp[,whichint]),,drop=F]
    diffmz=diff(xtmp[,whichmz])
    while(any(diffmz<deltaMZ) & nrow(xtmp)>1){
      (imerg=which.min(diffmz))
      if(!keepMax){
        xtmp[imerg+1,whichmz]=weighted.mean(xtmp[imerg+(0:1),whichmz],xtmp[imerg+(0:1),whichint])
        xtmp[imerg+1,whichint]=sum(xtmp[imerg+(0:1),whichint])
      } else imerg=imerg+which.min(xtmp[imerg+(0:1),whichint])-1
      xtmp=xtmp[-imerg,]
      xtmp=xtmp[order(xtmp[,whichmz]),,drop=F]
      diffmz=diff(xtmp[,whichmz])
    }
    xtmp
  }
  
  toconv=c("E0"="low","E1"="mid","E2"="high")
  tmp=t(sapply(strsplit(scan(ifile,what="raw",sep="\n")," "),function(x) x[1:3]))
  if(tmp[1,1]=="ERROR") return(list(itmp['Id'],NULL))
  tmp[grep("ener",tmp[,1]),3]=tmp[grep("ener",tmp[,1]),1]
  for(i in which(is.na(tmp[,3]))) tmp[i,3]=tmp[i-1,3]
  tmp=data.frame(tmp[-grep("ener",tmp[,1]),])
  tmp[,1]=round(as.numeric(tmp[,1]),5)
  tmp[,2]=round(as.numeric(tmp[,2]),5)
  tmp=tmp[which(tmp[,2]>0),]
  tmp[,3]=toconv[gsub("energy","E",tmp[,3])]

  cmb=.infct(tmp)
  cmb[,3]="cmb"
  cmb[,2]=round(100*cmb[,2]/sum(cmb[,2]),6)
  cmb[,1]=round(as.numeric(cmb[,1]),5)
  tmp=rbind(tmp,cmb)
  tmp=tmp[tmp[,2]>0.0001,]
  
  
  ############
  # header
  header2exp=c('BEGIN IONS',
               sprintf("PEPMASS=%.5f", itmp['PrecMass']),
               sprintf("CHARGE=%s", itmp['ChargeTop']),
               'SEQ=*..*','NAME=',
               sprintf("IONMODE=%s", itmp['Mod']),
               paste0("EXACTMASS=",itmp['Mass']),
               #paste0("SMILES=",itmp$Smiles),
               'SCANS=1')
  ############
  # loop
  llexp=list()
  for(ien in names(which(table(factor(tmp[,3],unique(tmp[,3])))>=minpk))){
    m=as.matrix(tmp[tmp[,3]==ien,1:2])
    toexp=c(header2exp,apply(m,1,function(x) sprintf("%.5f %.4f",x[1],x[2])),"END IONS")
    toexp[5]=paste0("NAME=",itmp['Id'],":",ien)
    #    if(length(llexp)>0) toexp=c(,toexp)
    llexp[[ien]]=toexp
  }
  invisible(list(itmp['Id'],llexp))
}


#' Multiple CFM file to MGF
#'
#' This function loads a
#'
#' @param mf2exp Data frame
#' @param what Path to the input file
#' @param fileout Path to the input file
#' @param root Root file
#' @return List of list
#' @export
ManyCfm2Mgf<-function(mf2exp,what="Neg",fileout=NULL,root='/media/david/ef81e301-125c-45e8-8e74-9ed5bbd3fd3b/CFM'){
  err=c()
  lfi=paste0(root,'/',what,'/',substr(mf2exp$Id,1,2),"/",gsub(";","_",mf2exp$Id),".log")
  
  mf2exp$Mod=ifelse(what=="Pos","Positive","Negative")
  mf2exp$ChargeTop=ifelse(what=="Pos","1+","1-")
  mf2exp$PrecMass=mf2exp$Mass+ifelse(what=="Pos",1,-1)*1.007227
  mf2exp$PrecMass[mf2exp$Charge==ifelse(what=="Pos",1,-1)]=mf2exp$Mass[mf2exp$Charge==ifelse(what=="Pos",1,-1)]
  l2exp=(mf2exp$Charge%in%c(ifelse(what=="Pos",1,-1),0) & file.exists(lfi))
  if(any(!l2exp)) err=c(err,as.vector(mf2exp$Id[which(!l2exp)]))
  allmgf=lapply(which(l2exp),function(i) Cfm2Mgf(lfi[i],itmp=as.vector(mf2exp[i,]),minpk = 3))
  err=c(err,as.vector(unlist(lapply(allmgf,function(x) if(is.null(x[[2]])) return(x[[1]])))))
  if(is.null(fileout)) return(list(allmgf))
  top=unlist(lapply(allmgf,function(x) x[[2]]),rec=F)
  if(length(top)>1) top[2:length(top)]=lapply(top[2:length(top)],function(x) c("",x))
  if(length(top)>0){
    cat(unlist(top),file=paste0(fileout,".mgf"),sep="\n")
    top=top[names(top)=="cmb"]
    if(length(top)>0){
      if(top[[1]][1]=="") top[[1]]=top[[1]][-1]
      # fileout=paste0(root,"/CombCFM/",what,"/",imf,"-cmb.mgf")
      cat(unlist(top),file=paste0(fileout,"-cmb.mgf"),sep="\n")
    }
  }
  
  return(err)
}