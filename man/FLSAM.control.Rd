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
  \item{stck}{An \code{\link[FLCore]{FLStock}} object
}
  \item{tun}{An \code{\link[FLCore]{FLIndices}} object
}
  \item{default}{Specifies how the control object should be populated. Default option is "full", which takes a reasonable first guess at a configuration. Any other value will populate all matrices with NA values.
}
}
\details{
An \code{\link{FLSAM}} assessment requires three elements - an \code{\link[FLCore]{FLStock}} object, an \code{\link[FLCore]{FLIndices}} object and an \code{\link{FLSAM.control}} object. The FLSAM.control object contains the configuration options for the model, and is therefore key to the process. However, there are a great many configuration options and therefore setting up the control object is a relatively complex process. The purpose of the FLSAM.control function is to create a template (possibly containing sensible first guesses at a configuration) object that can then be populated by the user according to their wishes.

The key concept with the configuration of a SAM model is that some parameters are bound together, so that the model fits just a single parameter, but uses it in multiple places e.g. by using a single parameter for the catchabilites of ages 5-9 in a survey. We refer to this concept as "binding" the parameters together. 

Bindings are expressed in a matrix format in the FLSAM.control object. The rows of the matrix correspond to the fleets, whilst the column correspond to the ages - it is important to remember that while all age-fleet combinations are represented here, not all are valid i.e. a fleet may not contain information on all ages. Bindings are specified with integers - the integers indicate which parameter should be used in the specific fleet-age context represented by that point in the matrix. Repeating a number means that all those fleet-age combinations with the same number will be fitted with just one parameter i.e. bound together. Parameters that are not to be fitted, or are turned off, are specified by setting the value to NA.

FLSAM.control can populate the template in two different ways, as indicated by the "default" argument. If default="full", a reasonable first guess at a configuration is supplied - this consists of linear catchability models for the surveys and a fixed selectivity pattern in the fishery. Alternatively, if "default" is set to any other value, all values are set to NA and therefore turned "off".
}
\value{
Returns an FLSAM.control object with the following slots
\item{name}{A text field giving the name of the object}
\item{desc}{A text field describing the object}
\item{range}{A named vector containg key information about range covered by the object \itemize{
\item min - the youngest age modelled
\item max - the oldest age modelled (can be a plus group)
\item plusgroup - the agegroup considered as a plusgroup (if any)
\item minyear - the first year modelled
\item maxyear - the last year modelled
\item minfbar - the younger limit of ages used in calcualting the average fishing mortality, fbar
\item maxfbar - the older limit of ages used in calcualting the average fishing mortality, fbar
}}
\item{fleets}{A named vector indicating the type of data that it represents. 0 = catch. 1=<I don't know! "con" something>. 2=Numbers at age survey. 3. SSB index}
\item{plus.group}{Logical specifying whether the oldest group should be considered as a plus group}
\item{states}{Binds the fishing mortality random-walks together}
\item{logN.vars}{Binds the variances of the numbers-at-age random-walks together}
\item{catchabilties}{Binds the catchabilities of the surveys together}
\item{power.law.exps}{Allows power-law catchability models to be implemented, and possibly bound together. In the case where a value is set to NA, a linear model will be employed for that fleet-age combination}
\item{f.vars}{Binds the variances of the fishing mortality random-walks together (as opposed to "states", which binds the value of the random walks i.e. to bind one random walk to multiple fishing-mortalties at age).}
\item{obs.vars}{Binds the variances (standard deviations) of the observations together}
\item{srr}{Stock recruitment relationships. Value options are 0= Random walk. 1= Ricker. 2=Beverton-Holt.}
}

\author{
Mark R. Payne and Niels Hintzen
}


\section{Warning }{It is important to realise that the default configuration returned by FLSAM.control is almost certaintly not the best configuration of your model - it is merely a first guess setup for your convenience. Similarly, there is no guarantee that this model will lead to a converged solution. It is the users responsibility to find the most appropriate ("best" configuration of the model.} 

\seealso{
\code{\link{FLSAM}}
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
