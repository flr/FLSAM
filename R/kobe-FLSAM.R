kobeFLSamFn=function(object,what=c("trks","smry","pts","ellipse")[1],prob=c(0.75,0.5,.25),nits=1000,bmsy=1,fmsy=1){
  
  require(ellipse)
  require(reshape)
  
  trmnl=c(max(seq(dim(vcov(object))[1])[dimnames(vcov(object))[[1]]=="ssb"]),
          max(seq(dim(vcov(object))[1])[dimnames(vcov(object))[[1]]=="fbar"]))
  
  currentCov=vcov(object)[trmnl,trmnl]
  
  ssb =params(object)[seq(dim(params(object))[1])[params(object)[,"name"]=="ssb"],2:3]
  fbar=params(object)[seq(dim(params(object))[1])[params(object)[,"name"]=="fbar"],2:3]
  
  trks=cbind(melt(data.frame(year=range(object)["minyear"]:range(object)["maxyear"],ssb),id="year"),
             melt(fbar)[,2])
  names(trks)[2:4]=c("quantity","stock","harvest")
  
  ellipse=as.data.frame(ellipse(currentCov,level=prob[1]))
  
  names(ellipse)=c("stock","harvest")
  ellipse=transform(ellipse, stock=stock+rev(ssb[,1])[1],harvest=harvest+rev(fbar[,1])[1])
  
  rnd=c(t(mvrnorm(nits,c(rev(ssb[,1])[1],rev(fbar[,1])[1]),currentCov)))
  rnd=as.data.frame(t(array(c(rnd),c(2,nits))))
  names(rnd)=c("stock","harvest")
  
  kb=object=llply(list(trks=trks,ellipse=ellipse,pts=rnd), function(x,bmsy,fmsy) transform(x, stock=stock/bmsy, harvest=harvest/fmsy), bmsy=bmsy, fmsy=fmsy)
  
  return(object)}

#### exported function
kobeFLSAM=function(object,what=c("trks","smry","pts","ellipse")[1],prob=c(0.75,0.5,.25),nits=1000,bmsy=1,fmsy=1){
  
  require(plyr)
  
  if (length(object)>1){
    res=kobeSamFn(object,what=what,prob=prob,nits=nits,bmsy=bmsy,fmsy=fmsy)
  } else {  
    res=mlply(object,   function(x,what=what,prob=prob,nits=nits,bmsy=bmsy,fmsy=fmsy)
      kobeSamFn(x,what=what,prob=prob,nits=nits,bmsy=bmsy,fmsy=fmsy),
              what=what,prob=prob,nits=nits,bmsy=bmsy,fmsy=fmsy)
    
    res=list(trks   =ldply(res, function(x) x$trks),
             pts    =ldply(res, function(x) x$pts),
             smry   =ldply(res, function(x) x$smry),
             ellipse=ldply(res, function(x) x$ellipse))
  }
  
  if (length(what)==1) return(res[[what]]) else return(res[what])}

#library(ggplot2)
#kb=kobeFLSamFn(sam)
#ggplot(kb$ellipse)+geom_path(aes(stock,harvest))+geom_path(aes(stock,harvest),data=subset(kb$trks,quantity=="value"))


