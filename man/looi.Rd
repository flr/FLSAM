\name{looi}
\alias{looi}
\alias{loo}
\alias{loi}
\title{
Run a sequence of assessments leaving one fleet (non-catch) in or out
}
\description{
Perform a sequence of assessments where each of the fleets are excluded or included from
the assessment run. It is also possible to run intermediate states where 2 or more
fleets are excluded or included from the assessment
}
\usage{
looi(stck, tun, ctrl, type="full")

loo(stck, tun, ctrl)
loi(stck, tun, ctrl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stck}{
      FLSAM object that specifies the stock characteristics
}
  \item{tun}{
      FLIndices object to use within the 'looi' assessment and out of which
      fleets will be included or dropped
}
  \item{ctrl}{
      FLSAM.control object that specifies the assessment settings
}
  \item{type}{Selects the type of LOOI run to perform. Options are
      \itemize{
      \item type: 'full' runs all combinations of leaving fleets in and out (may take a long time to run)
      \item type: 'LOI'  runs all combinations keeping only 1 fleet in
      \item type: 'LOO'  runs all combinations leaving only 1 fleet out.}
Default value is 'full'.}
}
\details{
The 'loo' and 'loi' methods call the 'looi' method with the appropriate type selected, and can be used as convenient shotcut}
\value{
An \code{\link{FLSAMs}} object is returned where the \code{names} slot contains the combination of fleets used within the assessments}
\author{
Niels T. Hintzen and Mark R. Payne
}
\seealso{
\code{\link{FLSAMs}}
}
\examples{
#- Load data
data(NSH)
data(NSH.sam)

#- Run a LOI
one.in <- loi(NSH,NSH.tun,NSH.ctrl)
plot(one.in)

#- Run a LOO
one.out <- loo(NSH,NSH.tun,NSH.ctrl)
plot(one.out)

}
