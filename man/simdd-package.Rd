\name{simdd-package}
\alias{simdd-package}
\alias{simdd}
\docType{package}
\title{
  simdd -- Simulation of various directional distributions
}
\description{
To simulate from various directional distributions,as 
  described in the paper "A new unified approach for the simulation 
  of a wide class of directional distributions"
}
\details{
\tabular{ll}{
Package: \tab simdd\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-02-11\cr
License: \tab GPL-2\cr
}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Maintainer: John Kent <j.t.kent@leeds.ac.uk>
}
\references{
For details, see the background paper "A new unified approach for the
simulation of a wide class of directional distributions".
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
%%%\keyword{ package }
\examples{
# simulate 20 random unit vectors from a Bingham distribution
# with concentration matrix diag(c(1,2,-3)).
X=rBingham(20,c(1,2,-3))
# check the efficiency of the simulation algorithm
attr(X,"summary")
}
