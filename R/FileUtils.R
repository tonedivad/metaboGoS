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
