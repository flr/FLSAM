\name{lr.test}
\alias{lr.test}
Log ratio test for FLSAM object
\title{
Perform a log-ratio test for all models specified
}
\description{
Compare model fit performance among the models specified based on log-likelihood ratio's.
Different in log-likelihood is calculated and degrees of freedom are given. P values
are specified assuming a chi-quared distribution.
}
\usage{
lr.test(FLSAM,...)
lr.test(FLSAMs)
}
\arguments{
  \item{FLSAM}{ FLSAM object}
  \item{FLSAMs}{FLSAMs object}
  \item{\dots}{
    In combination with FLSAM object: additional objects of the same type (need at least 1)
  }
}
\value{
A summary table of the test characteristics is returned}
\warnings{Log-likelihood tests can only be used for nested models (i.e. model parameterisation is different)
and won't work for different input data or different model types.
}
\author{
Code by Anders Nielsen, implemented by Niels T. Hintzen and Mark R. Payne
}
\seealso{
\code{\link{looi}}
}
\examples{
#- Load the data
data(NSH.sam)

#- Run lr ratio test on FLSAMs object
lr.ratio(NSH.sams)

#- Run lr ratio test on FLSAM objects
lr.ratio(NSH.sam,NSH.sam,NSH.sam)
}