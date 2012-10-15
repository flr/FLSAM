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
params(object,...)
coef(object,...)
coefficients(object,...)
}
\arguments{
  \item{object}{An \code{\link{FLSAM}} or \code{\link{FLSAMs}} object.}
  \item{\dots}{Other arguments}
}
\details{
These accessor functions provide direct access to the values returned by the SAM model. However, interpretation of these parameters can be challenging and requires some knowledge of the underlying workings of the SAM model. For more access in a more user-friendly format, see the \code{\link{accessors}} help.

The following parameter accessor functions are currently supported: \itemize{
\item params - returns all values estimated by the model, including derived parameters such as SSB, TSB, etc.
\item coef - returns the parameters fitted by the model, but not derived parameters such as SSB, TSB etc.
\item coefficients - an alias for coef
}
Accessor functions are provided for both \code{\link{FLSAM}} and \code{\link{FLSAMs}} objects.

Note that the naming of these parameters is somewhat obtuse. A brief description of these parameters follows:\itemize{ 
\item{logFpar}{: Log-transformed catchability of the numbers-at-age survey indices. The binding of these parameters are controlled in the 'catchabilities' slot of the FLSAM.control object.}
\item{logSdLogFsta}{: Log-transformed standard deviation of the random walk in the (log-transformed) fishering mortalities. Bindings are controlled in the 'f.vars' slot.}
\item{logSdLogN}{: Log-transformed standard deviation of the random walk in (log-transformed) N. Bindings are controlled in the logN.vars slot}
\item{logSdLogObs}{: Log-transformed observation variance for each data source. Bindings are controlled in the 'obs.vars' slot}
\item{logScaleSSB}{: Log-transformed catchability of the SSB indices}
\item{logSdSSB}{: Log-transformed observation variance of the SSB indices}
\item{rec_loga}{: Log-transformed 'a' parameter of the stock-recruitment model. Fitting of this parameter is controlled by the 'srr' slot}
\item{rec_logb}{: Log-transformed 'b' parameter of the stock-recruitment model. Fitting of this parameter is controlled by the 'srr' slot}
\item{rho}{: Correlation coefficient in the correlated random-walk in fishing mortality. Fitting of this parameter is controlled by the 'cor.F' slot}
\item{logPowSSB}{: Log-transformed power-law coefficient of the SSB catchability model.}
\item{U}{: The state variables. The order of the state variables is N(year 1),F(year 1, N(year 2), F(year 2), where N are the numbers-at-age state variables and F the fishing mortality state variables)}}

Not all parameters may be present in a fit, depending on the specific configuration of the FLSAM.control.
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
ssb(NSH.sam)
fbar(NSH.sam)
rec(NSH.sam)
#Power law parameters could be extracted as follows, but this will
#return an error in this case, as there are no power-law exponents in the model
#power.law.exps(NSH.sam) 

#And for FLSAMs
data(HERAS.sams)
ssb(HERAS.sams)
rec(HERAS.sams)
fbar(HERAS.sams)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
