

.MGmultiEICintgr<-function(eic,parDeco=list(psdrt=0.00446,span=7,
                                            minHeightMS1=100000,minNoiseMS1=5000),l2excl=NULL,doPlot=F){

  m=t(eic$y2)#/1000000
  m0=t(eic$y)#/1000000
  mmz=t(eic$mz)#/1000000
  
  span=parDeco$span
  minNoise=parDeco$minNoiseMS1
  
  nspan=floor(span/2)*2+1
  
  
bsllist=list(bsl=m,bslsc=m)

for(isid in colnames(m)){
  y=m[,isid]
  y=c(rep(y[1],nspan*5+1),y,rep(rev(y)[1],nspan*5+1))
  bsl=GRMeta:::.GRbslrf(1:length(y),y,NoXP = NULL)
  bsl$fit[bsl$fit<minNoise]=minNoise
  bslscore <- (y - bsl$fit)/max(bsl$sigma, 10^-3)
  bslscore[which(abs(bslscore) > 10)] = sign(bslscore[which(abs(bslscore) > 10)]) * 10
  bsllist$bsl[,isid]=bsl$fit[nspan*3+2+(1:nrow(m))]
  bsllist$bslsc[,isid]=bslscore[nspan*3+2+(1:nrow(m))]
}

# pks=.MGsimpleIntegr(lrt,y,noise.local =bsl$fit,snr.thresh = 2,span=floor(span/2)*2+1,minNoise = minNoise*1.01)
# pks$Sid=i



## get the master peak list
#tabres=.GRcompSeg(m,bsllist$bsl,list(nspan=nspan),5000)
apks=list()
for(isid in colnames(m)){
  re=.MGsimpleIntegr(x=1:nrow(m),y=m[,isid],noise.local =bsllist$bsl[,isid],snr.thresh = 2,span=nspan,minNoise = minNoise*1.01,v2alim = 0.8)
  if(nrow(re)==0) next
  re$Sid=isid
  apks[[isid]]=re
}
apks=do.call("rbind",apks)
names(apks)=gsub("pk\\.","tick.",names(apks))
apks$PkCl=NA

# llpk=.GRisover(apks$tick.left,apks$tick.right,retOne = T)
# llpk=lapply(llpk,function(x) .GRsplist(apks$tick.loc[x],x,d=5*nspan))
llpk=.GRsplist(apks$tick.loc,d=6.1*nspan)
if(any(sapply(llpk,is.list))) llpk=unlist(llpk,recursive = F)
for(i in 1:length(llpk)) apks$PkCl[llpk[[i]]]=i

l2k=as.numeric(names(which(tapply(apks$tick.int*2,apks$Pk,max)>parDeco$minHeightMS1 & tapply(apks$tick.snr,apks$Pk,max)>2)))
l2k=l2k[l2k%in%unique(apks$PkCl[!apks$Sid%in%l2excl])]
apks=apks[apks$Pk%in%l2k,]
apks$PkCl=as.numeric(factor(apks$PkCl,names(sort(tapply(apks$tick.loc,apks$PkCl,mean)))))
apks=apks[order(apks$Sid,apks$PkCl),]
apks$PkCl2=1


l2chk=names(which(tapply(apks$tick.loc,apks$PkCl,function(x) diff(range(x)))>nspan |
                    tapply(apks$Sid,apks$PkCl,function(x) any(duplicated(x)))))
i2chk=l2chk[1]

for(i2chk in l2chk){
  l2verif=which(apks$PkCl==i2chk)
  apks$PkCl2[l2verif]=.infctverif(apks[l2verif,],m=m,nspan = nspan,l2excl = l2excl,doplot=(doPlot>1))
}
newcl=paste(apks$PkCl,apks$PkCl2,sep=";")
apks$PkCl2=as.numeric(factor(newcl,names(sort(tapply(apks$tick.loc,newcl,mean,na.rm=T)))))

if(doPlot>0){
  cols=rep(brewer.pal(8,"Dark2"),length(unique(apks$Sid))/8)
  cols=cols[as.numeric(factor(apks$Sid,unique(apks$Sid)))]
  plot(apks$tick.loc,col=cols,pch=21,cex=.5)
abline(h=tapply(apks$tick.loc,apks$PkCl2,mean),lty=2,col="grey")
text(1:nrow(apks),apks$tick.loc,apks$PkCl2,col=cols)
}
#####################
## Form the final peak table
llpks=tapply(1:nrow(apks),paste0(apks$PkCl2,";;",apks$Sid),list)

sc2rt=as.numeric(rownames(m))*parDeco$psdrt
minpk=is.na(m)
matpks=list()
for(idx in names(llpks)){
  i=llpks[[idx]]
  lx=range(apks[i,1:3])
  lx=lx[1]:lx[2]
  lrt=sc2rt[lx]
  minpk[lx,apks$Sid[i]]=T
  v=m[lx,apks$Sid[i]]
  vb=bsllist$bsl[lx,apks$Sid[i]]
  v0=m0[lx,apks$Sid[i]]
  vmz=mmz[lx,apks$Sid[i]]
  lnna=which(!is.na(v0))
  lnna=lnna[order(abs(lnna-which.max(v)),-v0[lnna])]
  if(length(lnna)>nspan) lnna=lnna[1:nspan] 
  iap=lnna[which.max(v0[lnna])]
  
  matpks[[idx]]=data.frame(Sid=apks$Sid[i][1],Pk=apks$PkCl2[i][1],
                           rtmin=min(lx),
                           rtmax=max(lx),
                           rt=lx[which.max(v)],
                           rtap=lx[iap],
                           npts=sum(!is.na(v0)),
                           mz=round(vmz[iap],6),
                           mzmed=round(weightedMedian(vmz,v0,na.rm=T),6),
                           snr.sm=round(max(v/vb),3),
                           int.sm=round(max(v),3),
                           intap=v0[iap],
                           area.sm=round(unname(.GRgetArea(lrt,v)[1]-.GRgetArea(lrt,vb)[1]),3),
                           area=round(.GRgetArea(lrt[!is.na(v0)],v0[!is.na(v0)])[1],3))
}

matpks=do.call("rbind",matpks)
## fix missing???
lupk=which(tapply(matpks$intap,matpks$Pk,max)>=parDeco$minHeightMS1 &
             tapply(matpks$snr.sm,matpks$Pk,max)>=2 &
             tapply(matpks$npts,matpks$Pk,max)>(nspan/2))
matpks=matpks[matpks$Pk%in%names(lupk),]
lupk=sort(tapply(matpks$rt,matpks$Pk,median))
matpks$Pk=as.numeric(factor(matpks$Pk,names(lupk)))
matpks=matpks[order(matpks$Pk,matpks$Sid),]

########################## Fix missing
pksum=cbind(round(tapply(matpks$rt,matpks$Pk,quantile,0.5)),
            floor(tapply(matpks$rtmin,matpks$Pk,quantile,0.25)),
            ceiling(tapply(matpks$rtmax,matpks$Pk,quantile,0.75)))
lmiss=names(which(table(matpks$Sid)!=max(matpks$Pk)))
for(isid in lmiss){
  lin=unique(matpks$Pk)[!unique(matpks$Pk)%in%matpks$Pk[matpks$Sid==isid]]
  for(ipk in lin){
    lx=pksum[ipk,2]:pksum[ipk,3]
    v=m[lx,isid]
    iv=minpk[lx,isid]
    if(all(v[!iv]<(minNoise*2))) next
    print(c(ipk,isid))
  }
}

### Set the rt back 
for(i in names(matpks)[grep("^rt",names(matpks))]) matpks[,i]=round(sc2rt[matpks[,i]],3)

rownames(matpks)=NULL
invisible(matpks)
}




.infctverif<-function(ipks,m,nspan,l2excl=NULL,doplot=FALSE){
  lsc=range(ipks[,1:3])+c(-1,1)*nspan
  lsc[1]=max(1,lsc[1]);lsc[2]=min(nrow(m),lsc[2]);lsc=lsc[1]:lsc[2]
  
  lusamp=unique(ipks$Sid)
  lusamp=lusamp[!lusamp%in%l2excl]
  re=nmf(m[lsc,lusamp], 1, 'snmf/l',seed='ica')
  #tabres=.GRcompSeg(m[lsc,isamp],bsllist$bsl[lsc,isamp],list(nspan=nspan),5000)
  
  linspan=rev(seq(nspan,nspan*2+1,2))
  inspan=1
  mpks=.GRmsPeakSimple(lsc,loadings(re)[,1],span = linspan[inspan])
  while(nrow(mpks)==0 & inspan<length(linspan)){
    inspan=inspan+1
    mpks=.GRmsPeakSimple(lsc,loadings(re)[,1],span = linspan[inspan])
  }
  linmpks=lapply(1:nrow(mpks),function(x) mpks[x,'mass.left']:mpks[x,'mass.right'])
  
  mperc=do.call("rbind",lapply(1:nrow(ipks),function(i){
    lsci=ipks[i,2]:ipks[i,3]
    sapply(linmpks,function(y) sum(m[lsci[lsci%in%y],ipks$Sid[i]]))/sum(m[lsci,ipks$Sid[i]])
  }))
  v=apply(mperc,1,which.max)
  v[apply(mperc,1,max)<0.1]=NA
  if(!doplot) return(v)
  
  plot(ipks$tick.loc,col=brewer.pal(11,"Spectral")[as.numeric(cut(apply(mperc,1,max),seq(0,1,.1)))],
       pch=15+apply(mperc,1,which.max),ylim=range(lsc))
  text(1:nrow(ipks),ipks$tick.loc,apply(mperc,1,which.max))
  abline(h=tapply(ipks$tick.loc,ipks$PkCl,mean))
  abline(h=mpks$mass.loc,col=2)
  
  invisible(v)
}
