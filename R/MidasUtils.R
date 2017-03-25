######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
### Midas 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' Interface to Midas
#'
#' Interface to Midas
#'
#' @param obj dda metaboGoS object
#' @param precmzrange range prec m/z
#' @param lsp2exp  names of the MS/MS spectra
#' @param limint Minimum relative intensity
#' @param tempdir Directory containing the output
#' @param basenam result file name
#' @return List of tempdir/file locations
#' @export
dda2Midas<-function(obj,precmzrange=NULL,lsp2exp=NULL,limint=10^-5,
                       tempdir="./",basenam=paste0("tmp",round(runif(1)*1000))){

  # tempdir="./"
  # basenam=paste0("tmp",round(runif(1)*1000))
  # limint=10^-5
  ofile=paste(tempdir,'/',basenam,".FT2",sep="")
  l2use=c()
  if(!is.null(lsp2exp)) l2use=obj$MS2Infos$SpId[which(obj$MS2Infos$SpId%in%lsp2exp)]
  if(length(precmzrange)==2) l2use=sort(unique(c(l2use,obj$MS2Infos$SpId[which(obj$RoiInfos$mzPrec>=precmzrange[1] & obj$RoiInfos$mzPrec<=precmzrange[2])])))
  data2exp=obj$MS2Infos[obj$MS2Infos$SpId%in%l2use,]
  data2exp$Sc=obj$MS2Infos[data2exp$SpId,]$Sc
  
  cat("H\tExtractor\tMeAgain",sep="\n",file = ofile,append = F)
  cat("H\tm/z\tIntensity\tResolution\tBaseline\tNoise\tCharge\n",sep="",file = ofile,append = T)
  cat("H\tInstrument Model\tTitof\n",sep="",file = ofile,append = T)
  cat(" --> Exporting ",nrow(data2exp)," spectra to ",ofile,"\n",sep="")
  for(i in 1:nrow(data2exp)){
    cat(sprintf("S\t%d\t%d\t%.5f\n",data2exp$Sc[i],data2exp$Sc[i],data2exp$mzPrec[i]),sep="",file = ofile,append = T)
    cat(sprintf("Z\t%d\t%.5f\n",ifelse(grepl("[pP]os",data2exp$Polarity[i]),1,-1),data2exp$mzPrec[i]),sep="",file = ofile,append = T)
    cat(sprintf("I\t%s\n",i,data2exp$SpId[i]),sep="",file = ofile,append = T)
    ispec=obj$MS2Data[[data2exp$SpId[i]]]
    ispec=ispec[order(-ispec[,"y"]),c("mz","y")]
    l=max(min(sum((ispec[,"y"]/max(ispec[,"y"]))>=limint),100),5)
    ispec=ispec[1:l,]
    ispec=ispec[order(ispec[,"mz"]),]
    cat(sprintf("%.5f\t%.1f\t0\t0\t0\t0",ispec[,"mz"],ispec[,"y"]),sep="\n",file = ofile,append = T)
  }
  cat("Done!!\n")
  invisible(list(tempdir=tempdir,basenam=basenam,ofile=ofile))
}

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
### Generate configuration file
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' Configuration file for Midas
#'
#' Configuration file for Midas
#'
#' @param metLib metabolite library
#' @param polarity polarity
#' @param precTol  precursor tolerance
#' @param fragTol fragment tolerance
#' @param breakRing break Ring TRUE/FALSE
#' @param fragDepth result file name
#' @param nSlaves num of CPU
#' @return string to print out
#' @export
genMidasCfg<-function(metLib,polarity="pos",precTol=0.1,fragTol=0.01,breakRing=TRUE,fragDepth=3,nSlaves=1){

c('[Metabolite_Identification]',
# Filename of the metabolite database
sprintf('Metabolite_Database = %s',normalizePath(metLib)),
# Presumed polarity if not provided in the input file. Options: "positive" and "negative"
sprintf('Default_Polarity = %s',ifelse(grepl("[pP]os",polarity),"positive","negative")),
# Presumed charge state for unknown charge state
'Default_Charge_State = 1',
##### Mass_Accuracy #####
# Mass Windows to be open around parent ion. -->
# Examples: a center window: "0", a center window and an offset window: "-1,0", etc -->
'Parent_Mass_Windows = 0\n',
'Positive_Ion_Fragment_Mass_Windows = 0, 1, 2',
'Negative_Ion_Fragment_Mass_Windows = -2, -1, 0',
# Parent mass tolerance in Da
sprintf('Mass_Tolerance_Parent_Ion = %.4f',precTol),
# Fragment mass tolerance in Da
sprintf('Mass_Tolerance_Fragment_Ions = %.4f',fragTol),
# Break ring structures in a metabolite or not: "true" or "false"
sprintf('Break_rings = %s',ifelse(breakRing,"true","false")),
# Depth of fragmentation pathways. The valid range is from 1 to 5
sprintf("Fragmentation_Depth = %d\n",fragDepth),
# Number of processes the search should initiate. This should be less than the total number of CPU cores.
sprintf("Number_of_Processes = %d\n",nSlaves))

}


######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
### Midas back to R
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
parseAnnotMone<-function(rootdir=".",rootnam="tmp"){
  
  ifile1=paste(rootdir,"/",rootnam,".FT2",sep="")
  ifile2=paste(rootdir,"/",rootnam,".AFT2",sep="")
  
  # ### load the original spectra
  # tmp=scan(ifile1,what="raw",sep="\n")
  # 
  # lst=grep("^[A-Z]",tmp)
  # lst=sapply(.GRsplist(lst,lst,d=1.01),min)
  # if(length(lst)>1) len=c(lst[-1]-1,length(tmp)) else len=length(tmp)
  # orispec=list()
  # for(i in 1:length(lst)){
  #   ispec=tmp[lst[i]:len[i]]
  #   sp=t(sapply(strsplit(ispec[grep("^[0-9]",ispec)],"\t"),function(x) as.numeric(x[1:2])))
  #   iisp=ifelse(any(grepl("^I",ispec)),"^I","^S")
  #   iisp=strsplit(ispec[grep(iisp,ispec)],"\t")[[1]][2]
  #   orispec[[iisp]]=sp[order(sp[,1]),,drop=F]
  # }
  # 
  ### load the annotated spectra
  tmp=scan(ifile2,what="raw",sep="\n",quiet = T)
  vn1=strsplit(tmp[1],"\t")[[1]]
  vn11=vn1[vn1%in%c( "ScanNumber","PrecursorMZ","MaximumMZ" ,"Rank","ParentMassError", "Score" ,  "ExplInt","MassMol"  )]
  vn2=strsplit(tmp[2],"\t")[[1]]
  vn21=vn2[vn2%in%c( "m/z","Intensity","NormIntensity" , "m/zError"   )]
  lst=grep("^M",tmp)
  if(length(lst)==0) return(NULL)
  if(length(lst)>1) len=c(lst[-1]-1,length(tmp)) else len=length(tmp)
  
  ameb=asp=list()
  for(i in 1:length(lst)){
    ispec=tmp[lst[i]:len[i]]
    ameb[[i]]=strsplit(ispec[1],"\t")[[1]][-1]
    toadd=data.frame(do.call("rbind",strsplit(ispec[-1],"\t")))
    colnames(toadd)=vn2
    for(ij in vn21) toadd[,ij]=as.numeric(toadd[,ij])
    names(toadd)[1:3]=c("x","y","ynorm")
    toadd$Ann=paste(toadd$SMILES,toadd$ProtonOffset,toadd$FragmentationLevel,sep=";")
    toadd$Ann[toadd$Ann=="NA;NA;NA"]=NA
    asp[[i]]=toadd[,c("x","y","ynorm","Ann","m/zError")]
  }
  ameb=data.frame(do.call("rbind",ameb))
  names(ameb)=vn1
  for(i in vn11) ameb[,i]=as.numeric(ameb[,i])
  rownames(ameb)=names(asp)=ameb$SpecId
  lso=order(-ameb$Score)
  allres=list(SpecAnnot=asp[lso],Matches=ameb[lso,])
  return(allres)
  
}
