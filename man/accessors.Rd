\name{accessors}
\alias{accessors}
\alias{catchabilities}
\alias{obs.var}
\alias{power.law.exps}
\alias{ssb}
\alias{tsb}
\alias{n}
\alias{rec}
\alias{f}
\alias{fbar}
\alias{catch}
\title{Accessor functions for FLSAM
}
\description{
Returns results from the FLSAM asessment model 
}
\usage{
catchabilities(object)
power.law.exps(object)

obs.var(object)

ssb(object)
tsb(object)

n(object)
rec(object)

f(object)
fbar(object)

catch(object)
}
\arguments{
  \item{object}{An \code{\link{FLSAM}} or \code{\link{FLSAMs}} :object.
}
}
\details{
Extracting many of the fitted values directly from the parameter slot of an \code{\link{FLSAM}} object can lead to difficulties in interpretation, due to the binding of parameters (i.e. fitting one parameter to multiple time series) and the fitting of parameters in log-space. These accessor functions provide easy access to these key results, correcting for these issues and providing both estimates of the parameter concerned and the associated confidence intervals.

The following accessor functions are currently supported: \itemize{
\item catchabilities - linear proportionality parameter relating the observations (e.g. survey indices) to the modelled values (e.g. numbers at age, SSB) 
\item power.law.exps  in cases where the relationship between the observations and the modelled values is modelled as being a power-law, this function returns the associated coefficients 
\item obs.var - the variances estimated internally within the model that are associated with each individual time series. These variances are mathematically analogous to weighting factors 
\item ssb - the spawning stock biomass estimated internally within the model 
\item tsb - the total stock biomass estimated internally within the model 
\item n - the estimate number of individuals at a given age in a given year 
\item rec - the estimated recruitment of fish at the youngest modelled age 
\item f - the fishing mortality on a given age in a given year 
\item fbar - the mean fishing mortality in a given year 
\item catch - the (modelled) total catch of fish (in weight) in a given year 
}
Accessor functions are provided for both \code{\link{FLSAM}} and \code{\link{FLSAMs}} objects.
}
\value{
The accessor functions return a \code{\link{data.frame}} whoose structure is dependent on the function called. Generally the data.frame can be considered to have two parts. 

The first part is the meta-data part, which gives context-specific information about the contents of each row. Possible columns returned  are \itemize{
\item name - In cases where an FLSAMs object is supplied, the name of the corresponding object is stored in this column. There is no requirement that FLSAMs object have unique names - however, extracting and processing data can be challenging if they do not!
\item fleet - The name of the fleet that the parameter corresponds to
\item year - The year that the parameter corresponds to
\item age - The age that the parameter corresponds to. In the case of biomass (or other non-age-based observations), age is NA.
}

The second part of the data-frame is common and constant to all accessor functions, and contains the actual parameters and their uncertainties \itemize{ 
\item value - The value of the parameter
\item CV - The coefficient of variation for the parameter estimate i.e. standard deviation of the log-transformed parameter
\item ubnd - The upper bound (95\% confidence limit) for the parameter estimate
\item lbnd - The lower bound (95\% confidence limit) for the parameter estimate
}}

\author{
Mark R. Payne
}

\seealso{
\code{\link{FLSAM}}, \code{\link{params}}, \code{\link{coef}}, \code{\link{coefficients}}
}
\examples{
#Load data
library(FLSAM)
data(NSH.sam)
#Extract parameter
catchabilities(NSH.sam)
obs.var(NSH.sam)
power.law.exps(NSH.sam) #Returns error, because no power-law exponents in model
ssb(NSH.sam)
fbar(NSH.sam)
rec(NSH.sam)
#And for FLSAMs
data(IBTS.sams)
ssb(IBTS.sams)
rec(IBTS.sams)
fbar(IBTS.sams)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
