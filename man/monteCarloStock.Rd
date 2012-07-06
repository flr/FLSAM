\name{monteCarloStock}
\alias{monteCarloStock}
\title{
Generate new starting conditions of numbers and f-at-age.
}
\description{
Method to sample from the variance-co-variance matrix and based on a multi-variante
normal distribution generate new starting conditions of numbers-at-age and f-at-age.}
\usage{
monteCarloStock(FLStock,FLSAM,realisations,run.dir=tempdir())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FLStock}{
An FLStock object
}
  \item{FLSAM}{
An FLSAM object
}
  \item{realisations}{
Number of newly generated starting conditions.
}
  \item{run.dir}{Path to which the temporary files are to be written to or read from. Default is the R temporary directory provided by tempdir()
}

}
\details{
Large number of realisations might take longer to create.
}
\value{
An FLStock with iter=number of realisations is returned}
\references{
}
\author{
N.T. Hintzen
}
\seealso{
\code{\link{FLSAM}}
}
\examples{
data(NSH)
data(NSH.sam)

summary(NSH)#iter = 1
newStock <- monteCarloStock(NSH,NSH.sam,100)
summary(newStock)#iter = 100
}