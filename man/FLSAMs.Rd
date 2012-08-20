\name{FLSAMs-class}
\Rdversion{1.1}
\docType{class}
\alias{FLSAMs-class}
\alias{FLSAMs}

\title{Class \code{"FLSAMs"}}
\description{
  \code{FLSAMs} is a class that stores multiple \code{FLSAM} objects in a single container, allowing quick and easy processing and plotting. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{FLSAMs(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{The data \code{"list"} ~~ }
    \item{\code{names}:}{Names of the elements \code{"character"} ~~ }
    \item{\code{desc}:}{Description of the object \code{"character"} ~~ }
    \item{\code{lock}:}{Lock mechanism. If turned on. the length of the list cannot be modified by adding or removing elements \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{FLlst}"}, list vector
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "FLSAMs", e2 = "FLStock")}: ... }
    \item{+}{\code{signature(e1 = "FLStock", e2 = "FLSAMs")}: ... }
    \item{catchabilities}{\code{signature(object = "FLSAMs")}: ... }
    \item{coefficients}{\code{signature(object = "FLSAMs")}: ... }
    \item{coef}{\code{signature(object = "FLSAMs")}: ... }
    \item{f}{\code{signature(object = "FLSAMs")}: ... }
    \item{FLSAMs}{\code{signature(object = "FLSAMs")}: ... }
    \item{lr.test}{\code{signature(object = "FLSAMs")}: ... }
    \item{obs.var}{\code{signature(object = "FLSAMs")}: ... }
    \item{power.law.exps}{\code{signature(object = "FLSAMs")}: ... }
	 }
}
\author{
Mark R. Payne and Niels Hintzen
}
\seealso{
\code{\link{accessors}}, \code{\link{HERAS.sams}}
}
\examples{
data(NSH)
data(NSH.sam)

#Modify HERAS catchability bindings systematically
#Note that this is the same code used to create the
#HERAS.sams data object
ctrl.list <- list(all.free=3:10,
                  four.lvls=c(3,3,4,4,5,5,6,6),
                  three.lvls=c(3,3,4,4,5,5,5,5),
                  two.lvls=c(3,3,3,4,4,4,4,4))
ctrls  <- lapply(ctrl.list,function(x) {
            tmp.ctrl <- NSH.ctrl
            tmp.ctrl@catchabilities["HERAS",ac(1:8)] <- x
            return(tmp.ctrl)})
sam.list <- lapply(ctrls,FLSAM,stck=NSH,tun=NSH.tun)                  
HERAS.sams <- FLSAMs(sam.list)

#Plot results
plot(HERAS.sams)

#Plot catchability coefficients for different levels
qs.dat <- catchabilities(HERAS.sams)
xyplot(value ~ age,data=qs.dat,groups=name,type="b",
    subset=fleet=="HERAS",auto.key=TRUE)  

#Update stock object
NSH.HERAS <- NSH + HERAS.sams
plot(NSH.HERAS)

#Compare models with AIC
HERAS.aics <-AIC(HERAS.sams)
plot(HERAS.aics,xaxt="n",xlab="Model",ylab="AIC")
axis(1,at=seq(HERAS.aics),labels=names(HERAS.aics)) #three lvls looks best
}
\keyword{classes}
