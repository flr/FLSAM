\name{AIC}
\alias{AIC}
\title{Akaike's information criterion for FLSAM
}
\description{
Returns Akaike's information criterion for a set of FLSAM objects
}
\usage{
AIC(object, ..., k = 2)
}
\arguments{
  \item{object}{An \code{\link{FLSAM}} or \code{\link{FLSAMs}} object
}
  \item{\dots}{Further \code{\link{FLSAM}}objects
}
  \item{k}{the _penalty_ per parameter to be used; the default
          'k = 2' is the classical AIC.
}
}
\details{
The AIC criteria can be used as a basis for comparing model fits in situations where a likelihood ratio test is not appropriate (i.e. the models are not nested). For more information see the help for \code{\link[stats]{AIC}} in the stats package.
}
\value{
Numeric value
}
\author{
Based on code by Anders Nielsen in SAM. Adapted to FLSAM by Mark R. Payne.
}

\seealso{
\code{\link{lr.test}}, \code{\link[stats]{AIC}}
}
\examples{
#Load assessment
library(FLSAM)
data(NSH.sam)
#Extract AIC
AIC(NSH.sam)
#For more an FLSAMs
data(IBTS.sams)
plot(AIC(IBTS.sams),pch=16)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

