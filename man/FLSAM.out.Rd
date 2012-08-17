\name{FLSAM.out}
\alias{FLSAM.out}
\title{
Create FLSAM formatted output file
}
\description{
Document the results of the FLSAM assessment including all input values, model settings and
output in formatted tables}
\usage{
FLSAM.out(FLStock, FLIndices, FLSAM, format = "TABLE \%i.")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FLStock}{
An FLStock object
}
  \item{FLIndices}{
An FLIndices object
}
  \item{FLSAM}{
An FLSAM assessment model output object
}
  \item{format}{
The Title to be specified with each table.
}
}
\value{
The output is a table that can be written to disk and pasted directly into
any text editor.
}
\author{
N.T. Hintzen
}
\examples{
data(NSH)
data(NSH.sam)

sam.out <- FLSAM.out(NSH,NSH.tun,NSH.sam,format="Table \%i")

#-Write to file, set margins
old.opt           <- options("width","scipen")
options("width"=80,"scipen"=1000)

write(sam.out,file=paste("C:/","sam.out",sep="."))
#-Reset margins
options("width"=old.opt$width,"scipen"=old.opt$scipen)

}
