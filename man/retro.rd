\name{retro}
\alias{retro}
\title{Retrospective stock assessment for FLSAM
}
\description{
Performs a retrospective stock assessment for the desired years using
the FLSAM stock assessment method
}
\usage{
retro(stock,indices,control,retro,year.range,return.FLStocks)
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
  \item{return.FLStocks}{ Logical. Defaul value = TRUE. If FALSE a list
                          containing FLSAM objects will be returned
}
}
\details{
The argument 'year.range' is a numeric vector of years for which the assessment is to be performed.
If this is not specified the integer ”retro” is used (default value = 0). If retro = 0
the assessment is run for final year only. If retro = 1, the assessment is run for the
penultimate and final year and so on.

The results of the restrospective analysis can be plotted using the generic FLStocks plot() method.}
\value{
Returns an object of type FLStocks or list of FLSAM objects. Each component FLStock object contains
the result of each of the retrospective assessments.
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

plot(res)
}

