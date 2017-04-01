### merge close MZ
infct.mergeclosemz<-function(currrois,currstats,thrMZ=0.001,dmz=0.0005,dppm=2.5){
  l=which(currstats[,"mndups"]>1)
  l=which(currrois$GrpEIC%in%currstats$GrpEIC)
  ll2chk=tapply(l,currrois[l,"GrpEIC"],c)
  l2merge=list()
  for(lx in ll2chk){
    tmproi=currrois[lx,]
    l2m=infct.chkclosemz(tmproi,thrMZ=thrMZ,dmz=dmz,dppm=dppm)
    if(length(l2m)==0) next
    l2merge=c(l2merge,lapply(l2m,function(y) lx[y]))
  }
  if(length(l2merge)==0) return(NULL)
  narois=do.call("rbind",lapply(l2merge,function(x) infct.merge(currrois[x,])))
  narois=rbind(currrois[-unlist(l2merge),],narois)
  narois=narois[order(narois$GrpEIC,narois$Sid,narois$intensity),]
  narois$GrpEIC=as.numeric(factor(narois$GrpEIC))
  nstats=infct.compstats(narois,"GrpEIC")
  invisible(list(rois=narois,stats=nstats))
}

###

infct.mergemzrt<-function(currrois,currstats,typ="mzmed"){
  
  if(length(typ)==1) csplit<-.GRsplistMZ(currstats[,mzmed],dppm=2.5,dmz = 0.0005)
  if(length(typ)==2) csplit<-.GRisover(currstats[,typ[1]],currstats[,typ[2]],ret=T)
  
  n=sapply(csplit,length)
  l2split=which(n>1)
  nsplit=list()
  cid=max(currstats[,"GrpEIC"])
  for(idx in l2split){
    lx=which(currrois$GrpEIC%in%currstats[csplit[[idx]],"GrpEIC"])
    tmproi=currrois[lx,]
    v=infct.loopmzrt(tmproi)
    newv=unlist(lapply(1:length(v),function(x) rep(x,length(v[[x]]))))
    if(sum(diag(table(newv,tmproi$GrpEIC)))==length(lx)) next
    nsplit=c(nsplit,list(cbind(x=lx,oGrpEIC=tmproi$GrpEIC,nGrpEIC=newv+cid)))
    cid=cid+max(newv)
  }
  if(length(nsplit)==0) return(NULL)
  nsplit=do.call("rbind",nsplit)
  currrois[nsplit[,1],"GrpEIC"]=nsplit[,3]
  narois=currrois[order(round(currrois$mz50),currrois$GrpEIC,currrois$Sid,currrois$intensity),]
  narois$GrpEIC=as.numeric(factor(narois$GrpEIC))
  nstats=infct.compstats(narois,"GrpEIC")
  invisible(list(rois=narois,stats=nstats))
}

#### Split based on rt/mz

infct.splitmzrt<-function(crois,cstats,lurois,cid=max(cstats[,"GrpEIC"])){
  nsplit=list()
  for(ieic in lurois){
    lx=which(crois$GrpEIC==ieic)
    tmproi=crois[lx,]
    v=.GRisover2(tmproi[,"mz10"],tmproi[,"mz90"],tmproi[,"rtmin"],tmproi[,"rtmax"],ret=T,thr1 = 0.001,thr2=0)
    if(length(v)<2) next
    newv=do.call("rbind",lapply(1:length(v),function(x) cbind(x=lx[v[[x]]],oGrpEIC=tmproi$GrpEIC[v[[x]]],nGrpEIC=x+cid)))
    nsplit=c(nsplit,list(newv))
    cid=max(newv[,3])
  }
  if(length(nsplit)==0) return(NULL)
  nsplit=do.call("rbind",nsplit)
 # l2rm=nsplit[,1]
  nrois=crois[nsplit[,1],]
  nrois$GrpEIC=nsplit[,3]
  #nrois$oGrpEIC=nsplit[,2]
  nstats=infct.compstats(nrois)
  luold=unique(nsplit[,2])
  ncstats=rbind(cstats[!cstats$GrpEIC%in%luold,],nstats)
  ncrois=rbind(crois[!crois$GrpEIC%in%luold,],nrois)
  return(list(rois=ncrois,stats=ncstats))
}

#### Split based on density mz

infct.splitmzdens<-function(crois,cstats,lurois,cid=max(cstats[,"GrpEIC"]),bw=5){
  nsplit=list()
  for(ieic in lurois){
    lx=which(crois$GrpEIC==ieic)
    tmproi=crois[lx,]
    v=.GRsplitOne(tmproi$mz,bw = bw,ismass = T)$sp
    if(length(v)<2) next
    newv=do.call("rbind",lapply(1:length(v),function(x) cbind(x=lx[v[[x]]],oGrpEIC=tmproi$GrpEIC[v[[x]]],nGrpEIC=x+cid)))
    nsplit=c(nsplit,list(newv))
    cid=max(newv[,3])
  }
  if(length(nsplit)==0) return(NULL)
  nsplit=do.call("rbind",nsplit)
  # l2rm=nsplit[,1]
  nrois=crois[nsplit[,1],]
  nrois$GrpEIC=nsplit[,3]
  #nrois$oGrpEIC=nsplit[,2]
  nstats=infct.compstats(nrois)
  luold=unique(nsplit[,2])
  ncstats=rbind(cstats[!cstats$GrpEIC%in%luold,],nstats)
  ncrois=rbind(crois[!crois$GrpEIC%in%luold,],nrois)
  return(list(rois=ncrois,stats=ncstats))
}

#####  #####  #####  #####  #####  #####  #####  #####  

infct.merge<-function(idf,orderby='intensity',
                      lminvar=c("rtmin","mz10","mzmin"),
                      lmaxvar=c("rtmax","mz90","mzmax","coda"),
                      lmedvar=c("mz50")){
  idf=idf[order(idf[,orderby],decreasing = T),]
  for(i in lminvar) idf[1,i]=min(idf[,i])
  for(i in lmaxvar) idf[1,i]=max(idf[,i])
  for(i in lmedvar) idf[1,i]=median(idf[,i])
  idf[1,]
}

#####  #####  #####  #####  #####  #####  #####  #####  

infct.chkclosemz<-function(tmproi,sid="Sid",mzmed="mz50",mzmin="mz10",mzmax="mz90",thrMZ=0.001,dmz=0.0005,dppm=2.5){
  l2merge=list()
  tab=table(tmproi[,sid])
  # print(tab[tab>1])
  lsid=names(which(tab>1))
  for(isid in lsid){
    l=which(tmproi[,sid]==isid)
    itmproi=tmproi[tmproi[,sid]==isid,]
    ill=lapply(.GRsplistMZ(itmproi[l,mzmed],dppm = dppm,dmz = dmz),function(x) l[x])
    ill=ill[lapply(ill,length)>1]
    if(length(ill)==0) next
    ill=unlist(lapply(ill,function(x) .GRisover(tmproi[x,mzmin],tmproi[x,mzmax],retOne = T,thr = thrMZ)),rec=F)
    ill=lapply(ill,function(x) l[x])
    l2merge=c(l2merge,ill[lapply(ill,length)>1])
  }
  return(l2merge)
}

#####  #####  #####  #####  #####  #####  #####  #####  

infct.loopmzrt<-function(tmproi,mzmin="mz10",mzmax="mz90",rtmax="rtmax",rtmin="rtmin",thrMZ=0.001,thrRT=0,maxiter=6){
  v=list(1:nrow(tmproi))
  nv=1;doLoop=TRUE;iter=0
  while(doLoop & iter<=maxiter){
    iter=iter+1
    v=.GRisover(tmproi[,mzmin],tmproi[,mzmax],retOne = T,thr = thrMZ)
    v=unlist(lapply(v,function(y){
      z=.GRisover(tmproi[y,rtmin],tmproi[y,rtmax],retOne = T,thr = thrRT)
      lapply(z,function(i) y[i])
    }),rec=F,use=F)
    doLoop=length(v)!=nv
    nv=length(v)
  }
  return(v)
}

#####  #####  #####  #####  #####  #####  #####  #####  

infct.compstats<-function(arois,grp="GrpEIC",mz="mz",mzmed="mz50",mzmin="mz10",mzmax="mz90",rtmin="rtmin",rtmax="rtmax"){
  # grp="GrpEIC";mzmed="mz50";mzmin="mz10";mzmax="mz90";rtmin="rtmin";rtmax="rtmax"
  # currsplit=tapply(1:nrow(arois),arois[,grp],c)
  # currstats=t(sapply(currsplit,function(x){
  #   v=10^6*(range(arois$mz50[x])/median(arois$mz50[x])-1)
  #   vrt=range(arois$rtmin[x],arois$rtmax[x])
  #   c(median(arois$mz50[x]),diff(v),v,median(arois$rt[x]),vrt,length(unique(arois$Sid[x])),max(table(arois$Sid[x])),length(x))
  # }))
  
  currstats=data.frame(do.call("rbind",tapply(1:nrow(arois),arois[,grp],function(x){
    v0=10^6*(range(arois[x,mz])/median(arois[x,mz])-1)
    v=10^6*(range(arois[x,mzmed])/median(arois[x,mzmed])-1)
    vrt=range(arois[x,rtmin],arois[x,rtmax])
    c(arois[x[1],grp],
      median(arois[x,mz]),diff(v0),
      median(arois[x,mzmed]),diff(v),v,
      range(arois[x,mzmin],arois[x,mzmax]),
      median(arois$rt[x]),vrt,length(unique(arois$Sid[x])),max(table(arois$Sid[x])),length(x))
  })))
  names(currstats)=c(grp,mz,"dppm0",mzmed,"dppm","dppmmin","dppmax","mz10","mz90","rt",rtmin,rtmax,"nsid","mndups","n")
  return(currstats)
}


# infct.compstats2<-function(arois,currsplit,mzmed="mz50",mzmin="mz10",mzmax="mz90",rtmin="rtmin",rtmax="rtmax"){
#   currstats=t(sapply(currsplit,function(x){
#     v=10^6*(range(arois$mz50[x])/median(arois$mz50[x])-1)
#     vrt=range(arois$rtmin[x],arois$rtmax[x])
#     c(median(arois$mz50[x]),diff(v),v,median(arois$rt[x]),vrt,length(unique(arois$Sid[x])),max(table(arois$Sid[x])),length(x))
#   }))
#   colnames(currstats)=c("mz","dppm","dppmmin","dppmax","rt","rtmin","rtmax","nsid","mndups","n")
#   return(currstats)
# }
