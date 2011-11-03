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
  \item{object}{An FLSAM object
}
  \item{\dots}{Further FLSAM objects 
}
  \item{k}{the _penalty_ per parameter to be used; the default
          ‘k = 2’ is the classical AIC.
}
}
\details{
The AIC criteria can be used as a basis for comparing model fits in situations where a likelihood ratio test is not appropriate (i.e. the models are not nested). For more information see the help for AIC in the stats package i.e. ?stats::AIC 
}
\value{
Numeric value
}
\author{
Based on code by Anders Nielsen in SAM. Adapted to FLSAM by Mark Payne
}

\seealso{
"lr.test", "stats::AIC"
}
\examples{
#Load assessment
library(FLSAM)
data(NSH.sam)
#Extract AIC
AIC(NSH.sam)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
