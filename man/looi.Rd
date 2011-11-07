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
looi(stck, tun, ctrl, type)
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
  \item{\dots}{\itemize{
      \item type: 'full' runs all combinations of leaving fleets in and out (may take a long time to run)
      \item type: 'LOI'  runs all combinations keeping only 1 fleet in
      \item type: 'LOO'  runs all combinations leaving only 1 fleet out.}
}
}
\value{
An FLSAMs object is returned where in the names(FLSAMs) returnes the combination of fleets used
within the assessment}
\references{
}
\author{
Niels T. Hintzen and Mark R. Payne
}
\seealso{
\code{\link{lr.ratio}},\code{\link{AIC}}
}
\examples{
#- Load data
data(NSH.sams)

#- Run a LOI
looi(NSH,NSH.tun,NSH.ctrl,type="LOI")

#- Run a LOO
looi(NSH,NSH.tun,NSH.ctrl,type="LOO")

#- Run a full evaluation of all combinations
looi(NSH,NSH.tun,NSH.ctr)
}
