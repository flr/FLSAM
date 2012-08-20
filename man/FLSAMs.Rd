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
\code{\link{accessors}}
}
\examples{
data(IBTS.sams)
data(NSH)

#Plot results
plot(IBTS.sams)

#Update stock object
NSH.new <- NSH + IBTS.sams

#Compare model fits
lr.test(IBTS.sams)
AIC(IBTS.sams)

}
\keyword{classes}
