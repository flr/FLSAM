\name{lr.test}
\alias{lr.test}
\title{Likelihood Ratio Test for FLSAM
}
\description{
Likelihood ratio test to compare the fit of two nested FLSAM assessment models
}
\usage{
lr.test(...)
}
\arguments{
  \item{\dots}{
FLSAM assessment objects. Currently limited to a maximum of two objects
}
}
\details{
A likelihood ratio test is an appropriate way to compare two nested models. 
}
\value{
A table detailing the results of the likelihood ratio test
}
\section{Warning! }{
lr.test() assumes that the user knowns what they are doing! A likelihood ratio test is only appropriate when the models in question can be considered as "nested" ie one model can be expressed as a subset of the other model. The lr.test() does not (and cannot) check for this nesting. In situations where the models to be compared cannot be considered to be nested, the use of the Akaike Information Criterion (AIC) is more appropriate.
} 
\author{
Based on code for SAM by Anders Nielsen. Implemented here by Mark Payne
}

\seealso{
"AIC"
}
\examples{
#None
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
