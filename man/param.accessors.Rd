\name{parameter.accessors}
\alias{catchabilities}
\alias{obs.var}
\alias{power.law.exps}
\title{Parameter Accessor functions for FLSAM
}
\description{
Returns some of the key parameters fitted by the FLSAM asessment model 
}
\usage{
catchabilities(object)
power.law.exps(object)
obs.var(object)
}
\arguments{
  \item{object}{An \code{\link{FLSAM}} object.
}
}
\details{
The observation variances, the survey catachabilites and any power-law exponents are some of the key parameters fitted by the SAM model. However, extracting these values directly from the parameter slot of an \code{\link{FLSAM}} object can lead to difficulties in interpretation, due to the binding of parameters (i.e. fitting one parameter to multiple time series). These accessor functions return a data-frame that corrects for these issues.

Note that these functions only return the corresponding value for indices at age - biomass-related parameters are not returned. This will hopefully be remedied in the future.
}
\value{
A \code{\link{data.frame}} with the following headers
\item{fleet}{The name of the fleet that the parameter corresponds to}
\item{age}{The age that the parameter corresponds to. In the case of biomass (or other non-age-based observations), age is NA.}
\item{param.name}{The name of the parameter as it is stored in the params slot}
\item{value}{The value of the parameter, "as fitted"}
\item{std.dev}{The standard deviation of the parameter}
}
\author{
Mark R. Payne
}
\note{
Parameters that are returned are the "as-fitted" parameters - any log-transformation etc that may be required is left up to the user. 
}

\seealso{
\code{\link{FLSAM}}
}
\examples{
#Load data
library(FLSAM)
data(NSH.sam)
#Extract parameter
catchabilities(NSH.sam)
obs.var(NSH.sam)
power.law.exps(NSH.sam) #Returns error, because no power-law exponents in model
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
