\name{drop.from.control}
\alias{drop.from.control}
\title{Drop an element from an FLSAM.control object
}
\description{
Drops multiple surveys or age groups from an FLSAM.control object
}
\usage{
drop.from.control(object, fleets="missing",ages="missing")
}
\arguments{
  \item{object}{An \code{\link{FLSAM.control}} object
}
  \item{fleets}{A vector containing the names of the fleets to be dropped
}
  \item{ages}{A vector containing the names of the age groups to be dropped
}
}
\details{
This function can be used to modify an \code{\link{FLSAM.control}} object by removing one or more fleets, or one or more age groups from the object. Changes in the sequence of bindings that may arise due to the removal of a fleet or age group are corrected for using the \code{\link[FLSAM]{update}} function.
}
\value{
An FLSAM.control object
}
\author{
Mark R. Payne.
}

\seealso{
\code{\link[FLSAM]{update}}, \code{\link{FLSAM.control}}
}
\examples{
#Load assessment
library(FLSAM)
data(NSH.sam)
#NSH.ctrl before dropping
NSH.ctrl
#Drop MLAI index
drop.from.control(NSH.ctrl,fleets="MLAI")
#Drop ages 8 and 9
drop.from.control(NSH.ctrl,ages=8:9)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
