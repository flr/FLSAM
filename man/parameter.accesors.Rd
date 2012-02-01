\name{parameter.accessors}
\alias{coef}
\alias{coefficients}
\alias{params}
\title{Parameter accessor functions for FLSAM
}
\description{
Returns parameters fitted by the FLSAM asessment model 
}
\usage{
params(object)
coef(object)
coefficients(object)
}
\arguments{
  \item{object}{An \code{\link{FLSAM}} or \code{\link{FLSAMs}} :object.
}
}
\details{
These accessor functions provide direct access to the values returned by the SAM model. However, interpretation of these parameters can be challenging and requires some knowledge of the underlying workings of the SAM model. For more access in a more user-friendly format, see the \code{\link{acessors}} help.

The following parameter accessor functions are currently supported: \itemize{
\item params - returns all values estimated by the model, including derived parameters such as SSB, TSB, etc.
\item coef - returns the parameters fitted by the model, but not derived parameters such as SSB, TSB etc.
\item coefficients - an alias for coef
}
Accessor functions are provided for both \code{\link{FLSAM}} and \code{\link{FLSAMs}} objects.
}
\value{
The parameter accessor functions return a \code{\link{data.frame}} whoose length depends on the calling function. params() returns all values (both fitted and dervied parameters), whereas coef() and coefficients() only return the fitted parameters.

Columns returned are \itemize{
\item object.name - In cases where an FLSAMs object is supplied, the name of the corresponding object is stored in this column. There is no requirement that FLSAMs object have unique names - however, extracting and processing data can be challenging if they do not!
\item index - An indexing number produced internally by SAM
\item name - The name of the parameter
\item value - The value of the parameter
\item std.dev - The standard deviation of the fitted parameter
}}

\author{
Mark R. Payne
}

\seealso{
\code{\link{FLSAM}}, \code{\link{accessors}}
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
