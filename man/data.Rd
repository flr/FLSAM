\name{NSH}
\alias{NSH}
\alias{NSH.tun}
\alias{NSH.ctrl}
\alias{NSH.sam}
\docType{data}
\title{
North Sea Autumn Spawning Herring Stock Assessment
}
\description{
A set of FLR objects used to perform the 2011 North Sea Autumn Spawning Herring Stock Assessment 
}
\usage{data(NSH)
data(NSH.sam)}
\format{
The data is split into two separate classes. The first of these, returned by \code{data(NSH)} contains objects relating to the North Sea herring stock
\itemize{
\item NSH - \code{\link[FLCore]{FLStock}} object for the North Sea herring stock, containing data relating to the stock e.g. catch numbers at age
\item NSH.tun - \code{\link[FLCore]{FLIndices}} object containing the results of surveys on the North Sea herring stock. The four surveys are 
  \itemize{
    \item MLAI - The Multiplicative Larval Abundance Index. Used as a index of SSB
    \item IBTS0 - A index of late-larval abundance performed in February. Used as an index of age 0 fish
    \item IBTS - The International Bottom Trawl Survey, performed in February. Used as an index of numbers-at-age for age 1-5 fish
    \item HERAS - The Herring Acoustic survey, performed in July. Used as an index of numbers at age for age 1-9 fish.
}}
In addition the results of an assessment of this stock performed using FLSAM are also included
\itemize{
\item NSH.sam - An \code{\link{FLSAM}} object containg the results of the assessment
\item NSH.ctrl - An \code{\link{FLSAM.control}} object containg the configuration used to perform the assessment
}
}

\details{
}
\source{
ICES 2011. Report of the Working Group for the Assessment of Herring South of 62N (HAWG).
}
\examples{
#Load data
data(NSH)
data(NSH.sam)
#Stock summary plot
plot(NSH.sam)
#Update NSH FLStock object
NSH <- NSH + NSH.sam
plot(NSH)

}
\keyword{datasets}
