\name{FLSAM.control}
\alias{FLSAM.control}
\title{FLSAM.control
}
\description{
Creates an FLSAM.control that can be used to control an FLSAM assessment
}
\usage{
FLSAM.control(stck, tun, default="full")
}
\arguments{
  \item{stck}{An FLSAM object
}
  \item{tun}{An FLIndices object
}
  \item{default}{Specifies how the control object should be populated. Default option is "full", which takes a reasonable first guess at a configuration
}
}
\details{
An FLSAM assessment requires three elements - an FLStock object, an FLIndices object and a FLSAM.control object. The FLSAM.control object contains the configuration options for the model, and is therefore key to the process. However, there are a great many configuration options and therefore setting up the control object is a relatively complex process. The purpose of the FLSAM.control function is to create a template (possibly containing sensible first guesses at a configuration) object that can then be populated by the user according to their wishes.

The key concept with the configuration of a SAM model is that some parameters are bound together, so that the model fits just a single parameter, but uses it in multiple places e.g. by using a single parameter for the catchabilites of ages 5-9 in a survey. We refer to this concept as "binding" the parameters together. The FLSAM.control object contains the definition bindings - 
}
\value{
Returns a FLSAM.control object with the following slots

}

\author{
Mark Payne and Niels Hintzen
}


 \section{Warning }{It is important to realise that the default configuration returned by FLSAM.control is almost certaintly not the best configuration of your model - it is merely a first guess setup for your convenience. Similarly, there is no guarantee that this model will lead to a converged solution.} ~

\seealso{
FLSAM
}
\examples{
#Load data
library(FLSAM)
data(NSH)
#Create an object
ctrl <- FLSAM.control(NSH,NSH.tun)
ctrl
#Bind HERAS catchabilities on ages 5:9
ctrl@catchabilities["HERAS",as.character(5:9)] <- 11 
#Alternatively can also use update functionality
ctrl@catchabilities["HERAS",as.character(5:9)] <- 101
ctrl <- update(ctrl)
#See example of full control object
data(NSH.sam)
NSH.ctrl

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
