\name{FLSAM.out}
\alias{FLSAM.out}
\alias{sam.out}
\title{
Create FLSAM formatted output file
}
\description{
Document the results of the FLSAM assessment including all input values, model settings and
output in formatted tables}
\usage{
FLSAM.out(stck, tun, sam, format = "TABLE \%i.")

sam.out(stck, tun, sam, format = "TABLE \%i.")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stck}{
An FLStock object to be documented
}
  \item{tun}{
An FLIndices object to be documented
}
  \item{sam}{
An FLSAM assessment model output object
}
  \item{format}{
The Title to be specified with each table.
}
}
\details{
This function produces a set of tables that can be used to directly document an assessment. Included in the tables are all (relevant) input data, the details of the tuning indices, the model parameter configuration, the estimated parameters and all estimated values. The structure and approach is inspired by the "ica.out" documentation from the ICA package.

\code{sam.out} is simply an alias for FLSAM.out.
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

#-Set margins to make prettier output
old.opt           <- options("width","scipen")
options("width"=80,"scipen"=1000)

#Write the file to the working directory
write(sam.out,file="sam.out")
#-Reset margins
options("width"=old.opt$width,"scipen"=old.opt$scipen)

}
