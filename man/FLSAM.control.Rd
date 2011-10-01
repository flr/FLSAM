\name{FLSAM.control}
\alias{FLSAM.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{FLSAM.control
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FLSAM.control(stck, tun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stck}{
%%     ~~Describe \code{stck} here~~
}
  \item{tun}{
%%     ~~Describe \code{tun} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(stck,tun) {
  #Default constructor
  #Create object
  ctrl <- new("FLSAM.control")

  #Populate metadata
  ctrl@range <- stck@range[c("min","max","minfbar","maxfbar")]
  ctrl@plus.group <- stck@range["plusgroup"]==stck@range["max"]
  
  #Setup coupling structures
  default.coupling <- matrix(as.integer(NA),nrow=1+length(tun),ncol=dims(stck)$age,
                        dimnames=list(fleet=c("catch",names(tun)),age=dimnames(stck@catch.n)$age))
  ctrl@states           <- default.coupling
  ctrl@catchabilities   <- default.coupling
  ctrl@power.law.exps   <- default.coupling
  ctrl@f.vars           <- default.coupling
  ctrl@obs.vars         <- default.coupling

  #Other variables
  ctrl@logN.vars            <- default.coupling[1,]
  ctrl@srr <- as.integer(0)

  #Finished!
  return(ctrl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
