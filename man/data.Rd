\name{NSH}
\alias{NSH}
\alias{NSH.tun}
\alias{NSH.ctrl}
\alias{NSH.sam}
\alias{IBTS.sams}
\docType{data}
\title{
North Sea Autumn Spawning Herring Stock Assessment
}
\description{
A set of FLR objects used to perform the 2012 North Sea Autumn Spawning Herring Stock Assessment 
}
\usage{data(NSH)
data(NSH.sam)
data(IBTS.sams)}
\details{
The data is split into three separate groups. The first of these, returned by \code{data(NSH)} contains objects relating to the North Sea herring stock
\itemize{
\item NSH - \code{\link[FLCore]{FLStock}} object for the North Sea herring stock, containing data relating to the stock e.g. catch numbers at age
\item NSH.tun - \code{\link[FLCore]{FLIndices}} object containing the results of surveys on the North Sea herring stock. The four surveys are 
  \itemize{
    \item MLAI - The Multiplicative Larval Abundance Index. Used as a index of SSB
    \item IBTS0 - A index of late-larval abundance performed in February. Used as an index of age 0 fish
    \item IBTS - The International Bottom Trawl Survey, performed in February. Used as an index of numbers-at-age for age 1-5 fish
    \item HERAS - The Herring Acoustic survey, performed in July. Used as an index of numbers at age for age 1-9 fish.
}}
In addition the results of an assessment of this stock performed using FLSAM, and the associated configuration file, are returned by \code{data(NSH.sam)}
\itemize{
\item NSH.sam - An \code{\link{FLSAM}} object containg the results of the assessment
\item NSH.ctrl - An \code{\link{FLSAM.control}} object containg the configuration used to perform the assessment
}
The third group is an example of an \code{\link{FLSAMs}} object and is called \code{IBTS.sams}. It contains five \code{\link{FLSAM}} objects resulting from a scan across the IBTS indices, successively tying the catchabilities and observation variances for the oldest ages together i.e. starting with all ages free, then binding age 4+, then 3+, then 2+ etc...

Note that these are not the full North Sea herring assessment objects, but rather have been truncated to remove the data prior to 1992, so as to allow the model to converge quickly. The full assessment can be found in the Herring Working Group repository, \link{http://hawg.googlecode.com}.
}

\source{
ICES. 2012. Report of the Herring Assessment Working Group for the Area South of 62 N (HAWG), 13 - 22 March 2012, Copenhagen, Denmark. ICES CM 2012/ACOM:06. 835 pp.
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
#Example of FLSAMs object
data(IBTS.sams)
#Likelihood ratio test to choose best combination
lr.test(IBTS.sams)

}
\keyword{datasets}
