\name{obsvar.plot}
\alias{obsvar.plot}
\title{
Create diagnostic observation variance plot.
}
\description{
Plot the estimated observation variance of each data source in the assessment
as a sorted barplot.
}
\usage{
obsvar.plot(FLSAM)
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
A barplot is returned with fleets at the horizontal and observation variance at
the vertical axis.}
\references{
}
\author{
Mark Payne, Niels T. Hintzen
}

\seealso{
\code{\link{obscv.plot}},\code{\link{residual.diagnostics}}}
}
\examples{
data(NSH.sam)

obsvar.plot(NSH.sam)
}
