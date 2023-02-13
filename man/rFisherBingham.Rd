\name{Simulation methods}
%\name{rFisherBingham}
\alias{rFisherBingham}
\alias{rBingham}
\alias{rBingham.Grassmann}
\alias{rFisher.SO3}
\alias{rvMsin.torus}
\alias{rBessel}

\title{Simulate Fisher Bingham and related distributions}

\description{
  Simulate the Fisher Bingham distribution on the unit sphere,
  the matrix Bingham distribution on a Grassmann manifold,
  the matrix Fisher distribution on SO(3), and the bivariate
  von Mises sin model on the torus = product of two circles.
}
\usage{
rFisherBingham(nsim, mu = 0, Aplus = 0, q = dimset(mu, Aplus), mtop = 1000)
rBingham(nsim, Aplus, q = dimq(Aplus), mtop = 1000)
rBingham.Grassmann(nsim, Aplus = 0, q = dimq(Aplus), r = 1, mtop = 1000)
rFisher.SO3(nsim, F)
rvMsin.torus(nsim, k1, k2, alpha, mtop=1000)
rBessel(nsim, k1, k2, alpha, mtop=1000)
}
\arguments{
   \item{nsim}{the desired number of simulations
%%     ~~Describe \code{nsim} here~~
}
\item{mu}{a vector of length q.  If mu is entered as a scalar, it
  is interpreted as the vector (0,...,0,mu) with zeros except in the
  final coordinate.
%%     ~~Describe \code{mu} here~~
}

\item{Aplus}{a symmetric q by q matrix.  If Aplus is given as a scalar,
  it is interpreted as the q by q matrix of zeros.  If Aplus is given as
  a vector, it is interpreted as the matrix diag(Aplus).}

\item{q}{The simulation is done on the unit sphere in R^q.  Thus for the
  circle, q=2.  Note that q only needs to be given explictly if 
  both mu and Aplus are scalars.  In the default setting
  q=dimset(mu,Aplus), the internal function dimset determines q from the
  length of mu or the dimension of Aplus, if it can.}

\item{r}{For rBingham.Grassmann, r denotes the number of columns (1 <=r
<=q-1) in the matrix Bingham distribution}

\item{mtop}{An internal parameter, designed to prevent infinite loops.}

\item{F}{The 3 by 3 parameter matrix for the matrix Fisher
      distribution in the function rFisher.SO3.}
    
\item{k1, k2, alpha}{The two concentration parameters and the
      interaction parameter in  the function rvMsin.torus.}
}

\details{
   The Fisher Bingham distribution on the unit sphere in R^q has density
  proportional to

  exp(mu^T x + x^T Aplus x)

  where x is a unit vector in R^q, and mu (q-dimensional vector) and
  Aplus (q by q symmetric) are parameters.

  The matrix Bingham distribution on q by r matrices X whose columns are
  orthonormal, is given by the density proportional to 

  exp(trace(X^T Aplus X)).

  The Bingham distribution on the unit sphere in R^q can be simulated
  using (a) rBingham, (b) rFisherBingham with mu=0, and (c)
  rBingham.Grassmann with r=1.  Choice (a) will be fastest.

  The Fisher distribution can be simulated using rFisherBingham with
  Aplus=0.

  The matrix Fisher distribution on SO(3) has density proportional to 

  exp(trace(F^T X))

  where X is a 3 by 3 rotation matrix, and F is a 3 by 3 parameter
  matrix.

  The bivariate von Mises sine model on the torus has density
  proportional to

  exp(k1 cos(theta) + k2 cos(phi) + alpha sin(theta) sin(phi))

  for two angles theta, phi.  The Bessel density is obtained from the
  bivariate von Mises sine model as the marginal density of theta.}

\value{For rBingham and rFisherBingham, the output is an nsim by q
  matrix.  Each row is a simulated unit vector.

  For rBingham.Grassmann, the output is an nsim by q by r array.  For
  each value of the first index, the result is a simulated q by r matrix with
  orthonormal columns.

  For rFisher.SO3,  the output is an nsim by r by r array.  For
  each value of the first index, the result is a simulated  r by r matrix
  rotation matrix.

  For rvMsin.torus, the output is an nsim by 4 matrix, with each row
  containing the simulated value of (cos(theta), sin(theta), cos(phi),
  sin(phi)).

  For rBessel, the output is an nsim by 2 matrix containing the marginal
  simulated values of (cos(theta), sin(theta)).

  In all cases, the output has an attribute "summary", which is a vector
  of length 6 summarizing some details about the number of simulations
  needed in the acceptance/rejection algorithm.  The key element of this
  vector is called "efficiency", a number between 0 and 1, where 1 means
  that all the simulated values from the envelope distribution have been
  accepted.  }

\references{
  For details, see the background paper "A new unified approach for the
simulation of a wide class of directional distributions".

%% ~put references to the literature/web site here ~
}
\author{J T Kent
%%  ~~who you are~~
}
%%%\note{
%%  ~~further notes~~
%%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  The function rvm in the CRAN library "CircStats" provides a more
  efficient method to simulate from the von Mises distribution (i.e. the
  Fisher distribution with q=2) than rFisherBingham.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
X1=rBingham(10,c(1,2,-3))
attr(X1,"summary")
X2=rFisherBingham(10,c(3,2,1),c(1,2,-3))
X3=rBingham.Grassmann(10,c(1,2,-3),r=2)
X4=rFisher.SO3(10,matrix(c(1,2,3,4,2,1,0,1,2),3,3))
X5=rvMsin.torus(10,2.2,3.1,1.3)
X6=rBessel(10,2.2,3.1,1.3)
}
%%%\keyword{ ~kwd1 }
%%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
