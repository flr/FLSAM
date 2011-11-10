\name{retro}
\alias{retro}
\title{Retrospective stock assessment for FLSAM
}
\description{
Performs a retrospective stock assessment for the desired years using
the FLSAM stock assessment method
}
\usage{
retro(stock, indices, control, retro, year.range)
}
\arguments{
  \item{stock}{An \code{\link{FLSAM}} object
}
  \item{indices}{An  \code{\link{FLIndices}} object
}
  \item{control}{An \code{\link{FLSAM.control}} object
}
  \item{retro}{ An integer that specifies the number of retrospective years.
                Default value = 0. Only used if year.range is not specified.
}
  \item{year.range}{Numeric vector of years to perform the assessment
}
}
\details{
The argument 'year.range' is a numeric vector of years for which the assessment is to be performed.
If this is not specified the integer ”retro” is used (default value = 0). If retro = 0
the assessment is run for final year only. If retro = 1, the assessment is run for the
penultimate and final year and so on.}
\value{
Returns an FLSAMs object. Individual assessments that failed to converge are not included in this object, but a warning is issued in each case.
}
\author{
Based on code by Laurence Kell, modified by Niels T. Hintzen and Mark R. Payne
}

\seealso{
\code{\link{looi}}
}
\examples{
#Load assessment
library(FLSAM)
data(NSH)
data(NSH.sam)

#Run a retrospective analyses for 3 years
res <- retro(NSH,NSH.tun,NSH.ctrl,retro=3)

#Make an FLStocks object and plot the results
NSH.retro <- res + NSH
plot(NSH.retro)
}

