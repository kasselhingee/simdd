\name{simdd-package}
\alias{simdd-package}
\alias{simdd}
\docType{package}
\title{
  simdd -- Simulation of various directional distributions
}
\description{
  Simulate from various directional distributions using methods
  described in the paper 'A new unified approach for the simulation 
  of a wide class of directional distributions'.
}
\details{
\tabular{ll}{
Package: \tab simdd\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2023-02-13\cr
License: \tab GPL-2\cr
}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
John Kent <j.t.kent@leeds.ac.uk>
}
\references{
John T. Kent, Asaad M. Ganeiber & Kanti V. Mardia (2018) A New Unified Approach for the Simulation of a Wide Class of Directional Distributions, Journal of Computational and Graphical Statistics, 27:2, 291-301, DOI: 10.1080/10618600.2017.1390468 
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
