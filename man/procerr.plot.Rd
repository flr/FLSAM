\name{procerr.plot}
\alias{procerr.plot}
\title{
Produce diagnostics on process error in the SAM Model}
\description{
Generate plots that inform on process error in the SAM assessment model
}
\usage{
procerr.plot(stck, weight = "stock.wt", type = "n", rel = F)
}
\arguments{
  \item{stck}{
An FLStock object
}
  \item{weight}{
Option to choose between stock.wt and catch.wt
}
  \item{type}{
Option to choose between observation error related to N-at-age ('n'), mortality ('mort')  or total biomass ('tsb')
}
  \item{rel}{
Logical. Express error relative or absolute
}
}
\value{
A picture is returned}
\author{
Niels T. HIntzen (based on code from Einar Hjorleifsson}
}
\examples{
library(FLSAM)

data(NSH)

procerr.plot(NSH,type="n")
procerr.plot(NSH,type="n",rel=T)
procerr.plot(NSH,type="mort")
procerr.plot(NSH,type="mort",rel=T)
procerr.plot(NSH,type="tsb",weight="catch.wt",rel=F)
procerr.plot(NSH,type="tsb",weight="catch.wt",rel=T)
procerr.plot(NSH,type="tsb",weight="stock.wt",rel=F)
procerr.plot(NSH,type="tsb",weight="stock.wt",rel=T)
}
