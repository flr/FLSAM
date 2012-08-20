\name{obscv.plot}
\alias{obscv.plot}
\title{
Create diagnostic observation variance v.s. CV plot.
}
\description{
Plot the estimated observation variance of each data source in the assessment
agains the CV of the estimate.
}
\usage{
obscv.plot(FLSAM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FLSAM}{
An FLSAM assessment output object
}
}
\details{
This plot can be used as rough and ready way to gauge how "tightly" the SAM model has converged on the observation variances. In cases where the observation variances are poorly defined, the associated CV of the estimate will be large. In such cases, the model convergence should be examined in more detail.
}
\value{
A plot is returned with on the horizontal axis the observation variance of each
data source versus on the vertical axis the CV estimate associated with the data
source}
\author{
Mark Payne, Niels T. Hintzen
}
\seealso{
\code{\link{obsvar.plot}},\code{\link{residual.diagnostics}}
}
\examples{
data(NSH.sam)

obscv.plot(NSH.sam)
}
