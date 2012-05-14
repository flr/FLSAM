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
}
\value{
A plot is returned with on the horizontal axis the observation variance of each
data source versus on the vertical axis the CV estimate associated with the data
source}
\references{
}
\author{
Mark Payne, Niels T. Hintzen
}
\note{
Note that combinations of either high observation variance with low CV or low
observation variance with high CV are indications that the results need to be
carefully scrutinized.
}
\seealso{
\code{\link{obsvar.plot}},\code{\link{residual.diagnostics}}}
}
\examples{
data(NSH.sam)

obscv.plot(NSH.sam)
}
