\name{simdd-package}
\alias{simdd-package}
\alias{simdd}
\docType{package}
\title{
  simdd -- Simulation of Fisher Bingham and Related Directional Distributions
}
\description{
Simulation methods for the Fisher Bingham distribution on the unit sphere, the matrix Bingham distribution on a Grassmann manifold, the matrix Fisher distribution on SO(3), and the bivariate von Mises sine model on the torus.
 The methods use the first ever general purpose acceptance/rejection simulation algorithm for the Bingham distribution and are described fully by Kent, Ganeiber and Mardia (2018).
 These methods superseded earlier MCMC simulation methods and are more general than earlier simulation methods.
 The methods may be slower in specific situations where there are existing non-MCMC simulation methods (see Section 8 of Kent, Ganeiber and Mardia (2018) for further details).
}
\details{
\tabular{ll}{
Package: \tab simdd\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2023-02-15\cr
License: \tab GPL-2\cr
}
}
\author{
John Kent <j.t.kent@leeds.ac.uk>
}
\references{
John T. Kent, Asaad M. Ganeiber & Kanti V. Mardia (2018) A New Unified Approach for the Simulation of a Wide Class of Directional Distributions, Journal of Computational and Graphical Statistics, 27:2, 291-301, DOI: \doi{10.1080/10618600.2017.1390468} 
}
\examples{
# simulate 20 random unit vectors from a Bingham distribution
# with concentration matrix diag(c(1,2,-3)).
X=rBingham(20,c(1,2,-3))
# check the efficiency of the simulation algorithm
attr(X,"summary")
}
