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
obsvar.plot(sam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sam}{
An FLSAM assessment output object
}
}
\details{
Observation variances in the SAM model are analogous to the "weightings" employed in other model - in this case, they represent how noisy each of the individual tuning fleets are. Close examination of the values of these observation variances can given valuable insight into what the model is doing, and where it thinks the most information is.
}
\value{
A barplot is returned with fleets at the horizontal and observation variance at
the vertical axis.}
\author{
Mark Payne, Niels T. Hintzen
}

\seealso{
\code{\link{obscv.plot}},\code{\link{residual.diagnostics}}
}
\examples{
data(NSH.sam)

obsvar.plot(NSH.sam)
}
