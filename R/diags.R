#for (i in names(getSlots("FLSAM"))[-21]) slot(sm,i)=slot(NSH.sam,i)
setGeneric('diags', function(object, ...)
  standardGeneric('diags'))
 
setMethod("diags", signature(object="FLSAM"),
  function(object){
  
    res=residuals(object)
    names(res)[c(4:6)]=c("y","yHat","residual")
    res=ddply(res, .(fleet,age), transform, residualLag=c(NA,rev(rev(residual)[-1])))
    res=ddply(res, .(fleet,age), transform, 
                    qq=qqnorm(residual,plot=F)[c("x","y")])
    names(res)[8:9]=c("qqx","qqy")
    res=merge(res,as.data.frame(stock.n(object),drop=T))
    names(res)[10]="x"
  
    return(res)})