#' @name CFM
#' @aliases Cfm2Mgf
#' @aliases ManyCfm2Mgf
#' 
#' @title Convert CFM file to MGF
#' some stufff
#'
#' @param ifile Path to the input file
#' @param itmp Path to the input file
#' @param minpk Path to the input file
#' @return List of list
#' @rdname CFM
#' @export
Cfm2Mgf<-function(ifile,itmp,minpk=3,comment=NULL){
  
  .infct<-function(xtmp,deltaMZ=0.001,whichmz="X1",whichint="X2",keepMax=TRUE){
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
  tmp=t(sapply(strsplit(scan(ifile,what="raw",sep="\n",quiet = T)," "),function(x) x[1:3]))
  if(ncol(tmp)==0) return(list(itmp['Id'],NULL))
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
  cmb=cmb[which(100*cmb[,2]/sum(cmb[,2]) >0.01),,drop=F]
  cmb[,1]=round(as.numeric(cmb[,1]),5)
  tmp=rbind(tmp,cmb)
  tmp=tmp[tmp[,2]>0.0001,]
  
  
  ############
  # header
  header2exp=c('BEGIN IONS',
               sprintf("PEPMASS=%.5f", round(itmp['PrecMass'],5)),
               sprintf("CHARGE=%s", itmp['ChargeTop']),
               'SEQ=*..*','NAME=',
               sprintf("IONMODE=%s", itmp['Mod']),
               paste0("EXACTMASS=",round(itmp['Mass'],5)),
               paste0("MOLFORM=",itmp['Mod']),
  paste0("SMILES=",itmp$Smiles))
if(length(comment)) header2exp=c(header2exp,'COM=CFM')
  header2exp=c(header2exp,'SCANS=1')
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
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#' @rdname CFM
#'
#' @param mf2exp Data frame
#' @param what Path to the input file
#' @param fileout Path to the input file
#' @param root Root file
#' @return List of list
#' @export
ManyCfm2Mgf<-function(mf2exp,what="Neg",fileout=NULL,root,minpk=3,comment=NULL){
  err=c()
  lfi=paste0(root,'/',what,'/',substr(mf2exp$Id,1,2),"/",gsub(";","_",mf2exp$Id),".log")
  
  mf2exp$Mod=ifelse(what=="Pos","Positive","Negative")
  mf2exp$ChargeTop=ifelse(what=="Pos","1+","1-")
  if(!'PrecMass'%in%names(mf2exp)){
    mf2exp$PrecMass=mf2exp$Mass+ifelse(what=="Pos",1,-1)*1.007227*(mf2exp$Charge!=0)
  mf2exp$PrecMass[mf2exp$Charge==ifelse(what=="Pos",1,-1)]=mf2exp$Mass[mf2exp$Charge==ifelse(what=="Pos",1,-1)]
  }
#  l2exp=(mf2exp$Charge%in%c(ifelse(what=="Pos",1,-1),0) & file.exists(lfi))
  l2exp=(file.exists(lfi))
  if(any(!l2exp)) err=c(err,as.vector(mf2exp$Id[which(!l2exp)]))
  allmgf=lapply(which(l2exp),function(i) Cfm2Mgf(ifile = lfi[i],itmp=as.vector(mf2exp[i,]),minpk = minpk,comment=comment))
  err=c(err,as.vector(unlist(lapply(allmgf,function(x) if(is.null(x[[2]])) return(x[[1]])))))
  if(is.null(fileout)) return(list(allmgf))
  top=unlist(lapply(allmgf,function(x) x[[2]]),rec=F)
  if(length(top)>1) top[2:length(top)]=lapply(top[2:length(top)],function(x) c("",x))
  if(length(top)>0){
    cat(unlist(top[names(top)!="cmb"]),file=paste0(fileout,".mgf"),sep="\n")
    top=top[names(top)=="cmb"]
    if(length(top)>0){
      if(top[[1]][1]=="") top[[1]]=top[[1]][-1]
      # fileout=paste0(root,"/CombCFM/",what,"/",imf,"-cmb.mgf")
      cat(unlist(top),file=paste0(fileout,"Cmb.mgf"),sep="\n")
    }
  }
  
  return(err)
}