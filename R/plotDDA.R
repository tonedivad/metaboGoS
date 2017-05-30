
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' Plot DDA object
#' 
#' @param sfile sfile
#' @param blfile blfile
#' @param dmz dmz
#' @param dppm dppm
#' @param drt drt
#' @param winMS2 winMS2
#' @keywords internal
#' 
#' @export
plotOneDDA<-function(sfile,blfile=NULL,dmz=0.001,dppm=5,drt=1,winMS2=0.501){
  
  
if(!file.exists(sfile)) stop(paste0(sfile," does not exist"))
  
  obj=parseOneDDA(sfile,defIsolation = winMS2)
  xRaw=xcmsRaw(sfile)
xRawBl=NULL
  if(!is.null(blfile)){
    
    if(!file.exists(blfile)) print(paste0("Blank file ",blfile," does not exist"))
    xRawBl=xcmsRaw(blfile)
    
  }
  
#### #### #### ####
### merge close PrecIons
llppmrt=.GRsplitBoth2(obj$MS2Infos$PrecMZ,obj$MS2Infos$RT,dppm = dppm,dmz=dmz,typ="min",drt = drt*2.1)
llppmrt<-lapply(llppmrt,function(x) x[order(obj$MS2Infos$RT[x])])
ms2inf2=data.frame(do.call("rbind",lapply(llppmrt,function(x){
  mzr=range(obj$MS2Infos[x,]$PrecMZ)
  rtr=range(obj$MS2Infos[x,]$RT)
  c(mzmin=mzr[1],mzmax=mzr[2],rtmin=rtr[1],rtmax=rtr[2])
})))
ms2inf2=ms2inf2[order(ms2inf2$mzmin),]
cat("--> ",nrow(ms2inf2)," sets of precursor ions\n",sep="")

#### #### #### ####
llperc=seq(1,nrow(ms2inf2),length.out = 11)
for(iroi in 1:nrow(ms2inf2)){
  if(iroi%in%llperc) cat(iroi," ")
  mzr=range(ms2inf2[iroi,c("mzmin","mzmax")])
  rtr=range(ms2inf2[iroi,c("rtmin","rtmax")])
  l=which(obj$MS2Infos$PrecMZ>=ms2inf2[iroi,"mzmin"] & obj$MS2Infos$PrecMZ<=ms2inf2[iroi,"mzmax"] & 
            obj$MS2Infos$RT>=ms2inf2[iroi,"rtmin"] & obj$MS2Infos$RT<=ms2inf2[iroi,"rtmax"])
  if(length(l)==0) next
  ms2infos=obj$MS2Infos[l,]
  ms2data=obj$MS2Data[ms2infos$SpId]
  .MGinoneDDAplot(ms2infos,ms2data,xRaw,xRawBl,dmz=dmz,dppm=dppm,drt=drt,winMS2=winMS2)
}
invisible(ms2inf2)
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#' Internal function for plotting DDA object
#' 
#' @param ms2infos ms2infos
#' @param ms2data ms2data
#' @param xRaw xRaw
#' @param xRawBl xRawBl
#' @param dmz dmz
#' @param dppm dppm
#' @param drt drt
#' @param winMS2 winMS2
#' @keywords internal
#' 
#' @export
.MGinoneDDAplot<-function(ms2infos,ms2data,xRaw,xRawBl=NULL,dmz=0.001,dppm=5,drt=1,winMS2=0.501){
  
  # dmz=0.001;dppm=5;drt=1;winMS2=0.501
  
  # mzr=range(ms2inf2[iroi,c("mzmin","mzmax")])
  # rtr=range(ms2inf2[iroi,c("rtmin","rtmax")])
  
  rtr=range(ms2infos$RT)
  mzr=range(ms2infos$PrecMZ)
  
  imain=sprintf("Prec=[%.4f;%.4f] RT=[%.2f;%.2f]",mzr[1],mzr[2],rtr[1],rtr[2])
  
  mzrwide=mzr+c(-1,1)*(winMS2+dmz)
  mzrprec=mzr+c(-1,1)*dppm*mzr*10^-6
  
  imain2=sprintf("Isolation: %.4f - %.4f",mzrwide[1],mzrwide[2])
  imain3=sprintf("Precursor: %.4f - %.4f",mzrprec[1],mzrprec[2])
  
  
  rtrwide=rtr+c(-1,1)*drt
  rtrwide[1]=max(0,rtrwide[1])
  
  
  
  # rtrwide=c(max(lrt[1],rtrwide[1]),min(lrt[2],rtrwide[2]))
  
  .infct<-function(xr,mzrprec,mzrwide,rtrwide,sid=0){
    m=.GRrawMat(xr,mzrange =mzrwide,rtrange=rtrwide*60)
    m[,"rt"]=round(m[,"rt"],4)
    if(nrow(m)==0) return(list(NULL,NULL))
    m=cbind(m,In=(m[,"mz"]>=mzrprec[1] & m[,"mz"]<=mzrprec[2])*1,sid=sid)
    
    rts=tapply(m[,"rt"],m[,"rt"],unique)
    m2=cbind(rt=rts,y=tapply(m[,"y"],m[,"rt"],sum)[as.character(rts)],
             yin=tapply(m[m[,"In"]==1,"y"],m[m[,"In"]==1,"rt"],sum)[as.character(rts)],
             you=tapply(m[m[,"In"]!=1,"y"],m[m[,"In"]!=1,"rt"],sum)[as.character(rts)],sid=sid)
    rownames(m2)=rownames(m)=NULL
    return(list(m,m2))
  }
  
  
  re=.infct(xRaw,mzrprec,mzrwide,rtrwide,sid=1)
  m=rbind(re[[1]])
  mtic=rbind(re[[2]])
  
  if(!is.null(xRawBl)){
    re2=.infct(xRawBl,mzrprec,mzrwide,rtrwide,sid=2)
    m=rbind(re[[1]],re2[[1]])
    mtic=rbind(re[[2]],re2[[2]])
  }
  
  ######
  cols=c(brewer.pal(10,"Paired")[c(3:6)],"grey50","grey20") ## bl/samp out/in
  typ=c("p","p","p","b","b","b")
  cex=c(1,1,1,1,.5,.5,.5)
  pch=c(16,16,16,21,18,16)
  lwd=c(1,1,1,2,1,1.5)
  what=c("you","yin","you","yin","y","y")
  lll=list(which(mtic[,"sid"]==2 & !is.na(mtic[,"you"])),
           which(mtic[,"sid"]==2 & !is.na(mtic[,"yin"])),
           which(mtic[,"sid"]==1 & !is.na(mtic[,"you"])),
           which(mtic[,"sid"]==1 & !is.na(mtic[,"yin"])),
           which(mtic[,"sid"]==2 & !is.na(mtic[,"y"])),
           which(mtic[,"sid"]==1 & !is.na(mtic[,"y"])))
  
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  nf <- layout(matrix(c(1,2,3,4,4,4), 2, 3, byrow = TRUE), respect = TRUE)
  #  layout.show(nf)
  
  par(mar=c(4,6,2,.5))
  ###### Plot TICs
  xl=pretty(rtrwide)
  yl=pretty(sqrt(c(0,max(mtic[,"y"],na.rm=T))))
  yt=range(yl)/(median(yl)*10^-6)
  plot(range(rtrwide),range(yl),cex=0,ylim=range(yl)+c(0,diff(range(yl))/20),xlim=range(xl)+c(0,diff(range(xl))/7),xlab="RT",ylab="TIC [Sqrt]",axes=FALSE,main=imain)
  axis(1,at=xl,pos=0)
  axis(2,at=yl,las=2,pos=min(xl))
  l=which(sapply(lll,length)>0)
  l=l[!l%in%c(1,3)]
  for(i in l) points(mtic[lll[[i]],"rt"],sqrt(mtic[lll[[i]],what[i]]),col=cols[i],pch=pch[i],typ=typ[i],lwd=lwd[i])
  abline(v=ms2infos$RT)
  if(!is.null(xRawBl)) legend("topright",c("Prec. (S)","Win. (S)","Prec. (B)","Win. (B)"),col=cols[c(4,6,2,5)],pch=16,bty="n",ncol=1)
  if(is.null(xRawBl)) legend("topright",c("Prec. (S)","Win. (S)"),col=cols[c(4,6)],pch=16,bty="n",ncol=1)
  ########### plot MZs
  yl=pretty(range(mzrwide))
  
  plot(range(xl),range(yl),ylim=range(yl),cex=0,xlim=range(xl),xlab="RT",ylab="",axes=F,main=imain2)
  abline(h=median(mzrprec))
  abline(h=floor(yl/0.1)*0.1,lty=2,col="grey")
  axis(1,at=xl,pos=yl[1])
  axis(2,at=yl,las=2,pos=xl[1])
  points(m[m[,"sid"]==2,"rt"],m[m[,"sid"]==2,"mz"],col=cols[1],typ="p",pch=16,cex=.5)
  points(m[m[,"sid"]==1 & m[,"In"]==0,"rt"],m[m[,"sid"]==1 & m[,"In"]==0,"mz"],col=cols[3],typ="p",pch=16)
  points(m[m[,"sid"]==1 & m[,"In"]==1,"rt"],m[m[,"sid"]==1 & m[,"In"]==1,"mz"],col=cols[4],typ="p")
  # points(ms2infos$RT,ms2infos$PrecMZ,pch=17,col=4,cex=1.5)
  
  
  #abline(v=ms2infos[lprecs,c("RT")])
  
  ########### plot MZs
  yl=pretty(c(range(mzrprec),median(mzrprec)+c(-1,1)*median(mzrprec)*10^-6))
  yt=range(yl)/(median(yl)*10^-6)
  
  plot(range(xl),range(yl),ylim=range(yl),cex=0,xlim=range(xl),xlab="RT",ylab="",axes=FALSE,main=imain3)
  abline(h=(ceiling(yt[1]):floor(yt[2]))*10^-6*median(yl),lty=2,col="grey")
  axis(2,at=yl,las=2,pos=xl[1])
  axis(1,at=xl,pos=yl[1])
  points(m[m[,"sid"]==2,"rt"],m[m[,"sid"]==2,"mz"],col=cols[1],typ="p",pch=16,cex=.5)
  points(m[m[,"sid"]==1 & m[,"In"]==0,"rt"],m[m[,"sid"]==1 & m[,"In"]==0,"mz"],col=cols[3],typ="p",pch=16)
  points(m[m[,"sid"]==1 & m[,"In"]==1,"rt"],m[m[,"sid"]==1 & m[,"In"]==1,"mz"],col=cols[4],typ="p")
  spcols=rep(brewer.pal(9,"Set1"),ceiling(length(ms2data)/2))
  points(ms2infos$RT,ms2infos$PrecMZ,pch=17,col=spcols,cex=1.5)
  
  #  plot.new()
  ########### plot spectra
  for(i in 1:length(ms2data)) ms2data[[i]][,"sp"]=i
  sp=do.call("rbind",ms2data)
  sp=sp[order(sp[,2]),]
  sp=cbind(sp,y2=sqrt(sp[,"y"]))
  yl=pretty(c(0,sp[,"y2"]))
  yll=range(yl)
  yll[2]=yll[2]+diff(yll)/5
  mzstep=round(diff(range(sp[,"mz"]))/20)
  mzstep=max(1,ceiling(mzstep/5)*5)
  xl=round((floor(min(sp[,"mz"])/mzstep):ceiling(max(sp[,"mz"])/mzstep))*mzstep)
  plot(mean(sp[,"mz"]),mean(sp[,"y2"]),cex=0,ylim=yll,xlim=range(xl),axes=F,ylab="Sqrt [Int]",xlab="m/z",main=NULL)
  abline(v=median(mzr),lty=2,col="grey")
  axis(1,xl,pos=0);axis(2,yl,las=2,pos=min(xl))
  points(sp[,c("mz","y2")],col=spcols[sp[,"sp"]],pch=16)
  for(i in order(sp[,"y2"])) segments(sp[i,"mz"],0,sp[i,"mz"],sp[i,"y2"],col=spcols[sp[i,"sp"]])
  legend("top",as.character(round(ms2infos$RT,2)),bty="n",ncol=12,lwd=2,col=spcols)
  
  par(def.par)
}

  
  
  