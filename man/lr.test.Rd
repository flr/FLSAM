\name{lr.test}
\alias{lr.test}
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

\section{Warning!}{Log-likelihood tests can only be used for nested models (i.e. model parameterisation is different) and are not appropriate for different input data or different model types.
}
\author{
Code by Anders Nielsen, implemented by Niels T. Hintzen and Mark R. Payne
}
\seealso{
\code{\link{looi}}
}
\examples{
#- Load the data
data(IBTS.sams)

#- Run lr ratio test on example FLSAMs object
lr.test(IBTS.sams)

#- Run lr ratio test on FLSAM objects
lr.test(IBTS.sams[[1]],IBTS.sams[[2]])
}
