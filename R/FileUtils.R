######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' Read MGF file
#'
#' @param filein list of files 
#' @param namIdx append to spec nam
#' @param startIdx starting number for spectrum name
#' @return List of tempdir/file locations
#' @export

readMGFfile<-function(filein,namIdx='DBSpId',startIdx=1,minpk=1){
  
  infct<-function(mgf){
    desc.idx <- grep("=", mgf)
    desc <- mgf[desc.idx]
    spec <- mgf[-desc.idx]
    
    ms <- do.call(rbind, strsplit(spec, "[[:space:]]+"))
    mode(ms) <- "double"
    if (!length(ms)) ms <- matrix(numeric(), ncol = 2L)
    colnames(ms)=c("mz","y")
    r <- regexpr("=", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)), substring(desc,1L, r - 1L))
    return(list(desc,ms))
  }
  
  alldf=list()
  allsp=list()
  for(ifile in 1:length(filein)){
  if(!file.exists(filein[ifile])) next
  mgftmp <- scan(file = filein[ifile], what = "", sep = "\n", quote = "",allowEscapes = FALSE, quiet = TRUE)
  begin <- grep("BEGIN IONS", mgftmp) + 1L
  end <- grep("END IONS", mgftmp) - 1L
  if(!length(begin)) next
  
  res=lapply(1:length(begin),function(i) infct(mgftmp[begin[i]:end[i]]))
  names(res)=paste0(namIdx,(1:length(res))+startIdx-1)
  startIdx=startIdx+length(res)
  
  ludesc=unique(unlist(lapply(res,function(x) names(x[[1]]))))
  dfsp=do.call("rbind",lapply(names(res),function(x) c("SpecID"=x,res[[x]][[1]][match(ludesc,names(res[[x]][[1]]))])))
  dfsp=data.frame(dfsp)
  for(i in grep("MASS",dfsp)) dfsp[,i]=as.numeric(dfsp[,i])
  if(length(filein)>1) dfsp$FileName=ifile
  allsp[[ifile]]=lapply(res,function(x) x[[2]])
  alldf[[ifile]]=dfsp
}
  dfsp=GRMeta:::.joinDF(alldf,"SpecID")
  rownames(dfsp)=dfsp$SpecID
  lsp=unlist(allsp,rec=FALSE)
  lsp=lsp[!sapply(lsp,is.null)]
  lsp=lsp[sapply(lsp,nrow)>=minpk]
  dfsp=dfsp[which(dfsp$SpecID%in%intersect(dfsp$SpecID,names(lsp))),]
  lsp=lsp[rownames(dfsp)]
  
  return(list(DF=dfsp,MS=lsp))
}

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#' Export a list of files to MSP
#'
#' @param spec list of spectra 
#' @param plim min intensity in perc.
#' @param outfile file to be written
#' @param roundint rounding for intensity
#' @param roundmz rounding for m/z
#' @return List of tempdir/file locations
#' @export

export2MSP<-function(spec,outfile=NULL,plim=10^-5,roundint=0,roundmz=5,topSp=Inf){
  
  ## plim=0;roundint=0;roundmz=5
  specnam=names(spec)
  if(is.null(specnam)) specnam=paste0("SP",1:length(spec))
  
  idx=0
  pks=list()
  for(ix in 1:length(spec)){
    vmz=spec[[ix]][,"mz"]
    vint=spec[[ix]][,"y"]
    l2k=which(vint/max(vint)>plim)
    if(length(l2k)==0) next
    if(length(l2k)>topSp) l2k=l2k[rank(-vint[l2k])[1:topSp]]
    l2k=l2k[order(vmz[l2k])]
    idx=idx+1
    header2exp=c("",sprintf("Name:  %s",specnam[ix]),
                 sprintf("DB#:  %d", idx),
                 sprintf("Num Peaks:  %d", length(l2k)))
    
     spec2exp=sprintf(paste0("%.",roundmz,"f %.",roundint,"f; "),vmz[l2k],vint[l2k])
    spec2exp=sapply(split(spec2exp, ceiling(seq_along(spec2exp)/10)),paste,collapse = " ")
    pks[[ix]]=c(header2exp,spec2exp)
  }
  if(!is.null(outfile)) cat(unlist(pks),file=outfile,sep="\n")
  invisible(pks)
}



