#' Load DDA run from a mzXML file
#'
#' This function loads a DDA/DIA/AIF run from a mzXML file
#'
#' @param cefi Path to the input file
#' @param sid Sample id to given 
#' @param add2File Additional information to complement the file infos
#' @param defIsolation Default isolation window around precursor
#' @param save Save FALSE/TRUE
#' @param outfile Filename to save
#' @param stype Sample type
#' @param roundmz Rounding for the mass
#' @param verbose Print stuff out
#' @param useSid Add sid to the spectra names
#' @param rtlim Limit MS/MS to a certain range
##' @param MinIntBP Minimum MS/MS basepeak  height
#' @return List of list
#' @export
parseOneDDA<-function(cefi,sid=NULL,add2File=NULL,defIsolation=0.501,save=FALSE,outfile=NULL,stype=NA,roundmz=6,verbose=T,useSid=FALSE,rtlim=c(-Inf,Inf)){
  
  # sid=NULL;stype=NA;method="qn";roundmz=5;winsize=0.51;dppm=30;dmz=0.01;verbose=T;save=FALSE;outfile=NULL;useSid=FALSE;rtlim=c(.6,14.5);MinIntBP=1000
  
  ### only one CE !!!! -> to expand diff SP==1

  ##  ##  ##  ## File
  # completionTime=NA
  # cefiraw=paste0(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", cefi),c(".raw",".Raw",".RAW"))
  # cefiraw=cefiraw[file.exists(cefiraw)]
  # if(length(cefiraw)>0) completionTime=as.chron(file.info(cefiraw)[,"mtime"], format = c(dates = "Y-M-D", times = "h:m:s"))
  # 
  if(is.null(sid)) sid=gsub(".*/","",gsub("\\..*","",gsub(" ","",cefi)))
  if(nchar(sid)==0) sid=gsub(".*/","",gsub(" ","",cefi))
  if(is.na(stype)){
    if(grepl("^([blBLQCcSTDstda]+)_.*",sid)) stype=gsub("^([blBLQCcSTDstda]+)_.*","\\1",sid)
    if(grepl("^([blBLQCcSTDstda]+)_.*",cefi)) stype=gsub("^([blBLQCcSTDstda]+)_.*","\\1",cefi)
    if(grepl("/([blBLQCcSTDstda]+)_.*",cefi)) stype=gsub(".*/([blBLQCcSTDstda]+)_.*","\\1",cefi)
  }
  if(nchar(sid)==0){sid="unknwon";stype="unk"}
  if(sum(nchar(stype),na.rm=T)==0) stype="Unk"
  if(verbose) cat("Parsing the XML file '",cefi,"' for sample '",sid,"' of type '",stype,sep="")
  
  File=data.frame(Sid=sid,File=normalizePath(cefi))
  if(!is.null(add2File)) for(i in names(add2File)[!names(add2File)%in%names(File)]) File[,i]=add2File[i]
  rownames(File)=File$Sid
  
  ####################
  doc <- xmlRoot(xmlParse(cefi, useInternal = TRUE))
  tmp=xmlChildren(doc[[1]])

  softtype=NA
  convsoft=tmp[[which(names(tmp)=="dataProcessing")]]
  if(any(grepl("Proteo",xmlAttrs(xmlChildren(convsoft)$software)))) softtype="ProteoWizard"
  if(any(grepl("ReAdW",xmlAttrs(xmlChildren(convsoft)$software)))) softtype="ReAdW"
  
  if(is.na(softtype)) stop('Could not find the conversion software!')
  if(verbose)  cat(" converted by '",softtype,"'\n",sep="")
  
  .infct<-function(y){
    re=xmlAttrs(y)
    z=xmlChildren(y)
    if(!is.null(z$precursorMz)) re=c(re,xmlAttrs(z$precursorMz),precursorMz=xmlValue(z$precursorMz))
    re
  }
  
  lna=c('ScanMode',"num","msLevel",'retentionTime',
        'polarity','activationMethod','collisionEnergy','EnergyUnits','windowWideness',  "lowMz","highMz" ,
        'precursorCharge','precursorMz','precursorIntensity','basePeakIntensity','totIonCurrent','filterLine')
  ### ReadW
  if(softtype=="ReAdW") ascan=unlist(lapply(which(names(tmp)=="scan"),function(i){
    re=list(xmlAttrs(tmp[[i]]))
    x=xmlChildren(tmp[[i]])
    lsc=which(names(x)=="scan")
    if(length(lsc)>0) re=c(re,lapply(lsc,function(j) .infct(x[[j]])))
    re
  }),recursive = F)
  
  ### ProteoWizard
  if(softtype=='ProteoWizard') ascan=lapply(which(names(tmp)=="scan"),function(i) .infct(tmp[[i]]))
  
  ###
  ascan=lapply(ascan,function(x) x[lna])
  ascan=data.frame(do.call("rbind",ascan))
  names(ascan)=lna
  ascan$ScanMode=ifelse(ascan$msLevel!="1","ProductIon","Scan")
  ascan$polarity=ifelse(ascan$polarity=="-","Negative","Positive")
  ascan$retentionTime=round(floor(as.numeric(gsub("[A-Z]","",ascan$retentionTime))*1000/60)/1000,3)
  ascan$collisionEnergy=round(as.numeric(ascan$collisionEnergy),3)
  
  for(i in c("lowMz","highMz",'precursorMz')) 
    ascan[,i]=round(floor(as.numeric(ascan[,i])*10^(roundmz+1))/10^(roundmz+1),roundmz)
  for(i in c('num','precursorCharge','precursorIntensity','basePeakIntensity','totIonCurrent')) ascan[,i]=round(as.numeric(ascan[,i]),1)
  
  ascan=ascan[order(ascan$num),]
  ###
  
  ### ReadW
  if(softtype=="ReAdW") suppressWarnings(ascan$collisionEnergy<-as.numeric(gsub(".*@([a-z]+)([0-9\\.]+) .*","\\2",ascan$filterLine)))
  
  ###
  ascan$ScMS1=ascan$SpId=NA
  l=which(ascan$ScanMode=="Scan")
  for(i in 1:length(l)) ascan$ScMS1[l[i]]=i
  for(i in which(is.na(ascan$ScMS1))) ascan$ScMS1[i]=ascan$ScMS1[i-1]
  
  if(verbose) cat(' -- found ', sum(ascan$msLevel==1) ," MS1 and ", sum(ascan$msLevel==2) ," MS2\n",sep="")
  
  ### name windows
  toadd=ifelse(useSid,paste0("-",sid),"")
  ascan$SpId=sprintf("SP%.4f@%.3f%s",ascan$precursorMz,ascan$retentionTime,toadd)
  ascan$SpId[l]=sprintf("FS@%.3f%s",ascan$retentionTime[l],toadd)
  
  ### Load MS2
  tmp=mzR::openMSfile(cefi)
  msinfos=mzR:::header(tmp)
  if(!all(ascan$num==msinfos$seqNum)) stop('msinfos issues')
  
  sc2sp=tapply(ascan$ScMS1,ascan$num,c)
  apks=mzR::peaks(tmp)
  l2use=which(sapply(apks,nrow)>0)
  for(i in l2use) apks[[i]]=cbind(sc=i,mz=apks[[i]][,1],y=apks[[i]][,2],sp=sc2sp[i])
  
  #### Reduce based on RT/minInt MS2
  lsp=range(ascan$ScMS1[ascan$retentionTime>=rtlim[1] & ascan$retentionTime<=rtlim[2] & ascan$ScanMode=="Scan"])
  l2use=which(ascan$ScMS1%in%(lsp[1]:lsp[2]))
#  l2excl=which(ascan$BPInt<MinIntBP & ascan$ScanMode=="ProductIon")
  # l2excl=which(ascan$ScanMode=="ProductIon")
  # if(length(l2excl)>0) l2use=l2use[!l2use%in%l2excl]
  ascan=ascan[l2use,]
  
  ## MS2
  l=which(ascan$ScanMode=="ProductIon")
  infms2=ascan[l,]
  MS2Dat=tapply(1:nrow(infms2),infms2$SpId,function(x) do.call("rbind",apks[infms2$num[x]]))
  
  lv1=c("num" , "ScMS1","SpId", "msLevel","retentionTime","polarity","activationMethod","collisionEnergy",  'windowWideness',"lowMz","highMz", 
        "precursorCharge","precursorMz","precursorIntensity" ,"basePeakIntensity"  ,"totIonCurrent")
  
  lv2=c("Sc" , "ScMS1","SpId", "msLevel","RT","Polarity","activationMethod","CE",   'WinSize',"lowMz","highMz",
        "PrecCharge","PrecMZ","PrecInt" ,"BPInt"  ,"TIC")
  
  MS2Infos=infms2[,lv1]
  names(MS2Infos)=lv2
  MS2Infos$WinSize[which(is.na(MS2Infos$WinSize))]=defIsolation
  MS2Infos$WinSize=as.numeric(MS2Infos$WinSize)
  rownames(MS2Infos)=infms2$SpId

    ###### Finalise
  
  if(verbose) cat(' --> Final number of MS2 spectra: ', nrow(MS2Infos) ," between ",sprintf('%.3f-%.3f',min(MS2Infos$RT),max(MS2Infos$RT))," min.\n",sep="")
  
  ##  ##  ##  ## Scan2rt
  msinfos=mzR:::header(tmp)
  sc2rt=round(msinfos$retentionTime/60,5)
  Scan2rt=matrix(sc2rt,nrow=1,dimnames = list(sid,1:length(sc2rt)))
  
  ##  ##  ##  ##
  lso=order(MS2Infos$ScMS1,MS2Infos$RT,MS2Infos$PrecMZ)
  obj=list(File=File,MS2Infos=MS2Infos[lso,],MS2Data=MS2Dat[lso],Scan2rt=Scan2rt)
  class(obj) = append(class(obj), "ddaSet")
  if(save & is.null(outfile)) save(file=paste0(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", cefi),".rda"),obj)
  if(!is.null(outfile)) save(file=outfile,obj)
  invisible(obj)
  
}
