\name{monteCarloStock}
\alias{monteCarloStock}
\title{
Generate new starting conditions of numbers and f-at-age.
}
\description{
Method to sample from the variance-co-variance matrix and based on a multi-variante
normal distribution generate new starting conditions of numbers-at-age and f-at-age.}
\usage{
monteCarloStock(stck,sam,realisations,run.dir=tempdir())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stck}{An object of class \code{FLStock}
}
  \item{sam}{An object of class \code{FLSAM} 
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
\author{
N.T. Hintzen
}
\seealso{
\code{\link{FLSAM}}
}
\examples{
data(NSH)
data(NSH.sam)

#Stock with one iter
summary(NSH)#iter = 1
plot(NSH)
#Generate new stock with 100 iters
#newStock <- monteCarloStock(NSH,NSH.tun, NSH.sam,100, ncores=2L)
#Now does it look?
#summary(newStock)#iter = 100
#plot(newStock)  #Voila!

}
