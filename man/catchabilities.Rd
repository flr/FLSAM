\name{catchabilities}
\alias{catchabilities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessor functions for FLSAM
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
catchabilities(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object) 
{
    params <- subset(object@params, name == "logFpar")
    params$no <- 1:nrow(params)
    bindings <- as.data.frame(as.table(object@control@catchabilities), 
        responseName = "no")
    bindings <- subset(bindings, !is.na(bindings$no))
    res <- merge(params, bindings)
    res <- res[order(res$fleet, res$age), ]
    res$param <- res$name
    res$name <- ifelse(length(object@name) == 0, "<None>", object@name)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
