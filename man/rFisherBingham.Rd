\name{Simulate Fisher Bingham and related distributions}
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
rFisher.SO3(nsim, Fmat)
rvMsin.torus(nsim, k1, k2, alpha, mtop=1000)
rBessel(nsim, k1, k2, alpha, mtop=1000)
}
\arguments{
   \item{nsim}{the desired number of simulations
}
\item{mu}{a vector of length \code{q}.  If \code{mu} is entered as a scalar, it
  is interpreted as the vector \code{(0,...,0,mu)} with zeros except in the
  final coordinate.
}

\item{Aplus}{a symmetric \code{q} by \code{q} matrix.  If \code{Aplus} is given as a scalar,
  it is interpreted as the \code{q} by \code{q} matrix of zeros.  If \code{Aplus} is given as
  a vector, it is interpreted as the matrix \code{diag(Aplus)}.}

\item{q}{The simulation is done on the unit sphere in \eqn{R^q}.  Thus for the
  circle, \code{q=2}.  Note that \code{q} only needs to be given explicitly if 
  both \code{mu} and \code{Aplus} are scalars.  In the default setting
  the internal function \code{dimset}, \code{q=dimset(mu,Aplus)}, determines \code{q} from the
  length of \code{mu} or the dimension of \code{Aplus}, if it can.}

\item{r}{For \code{rBingham.Grassmann}, \code{r} denotes the number of columns (\code{1 <= r <= q-1}) in the matrix Bingham distribution}

\item{mtop}{An internal parameter, designed to prevent infinite loops.}

\item{Fmat}{The 3 by 3 parameter matrix for the matrix Fisher
      distribution in the function \code{rFisher.SO3}.}
    
\item{k1, k2, alpha}{The two concentration parameters and the
      interaction parameter in  the function \code{rvMsin.torus}.}
}

\details{
   The Fisher Bingham distribution on the unit sphere in \eqn{R^q} has density
  proportional to

  \deqn{\exp(\code{mu}^T x + x^T \code{Aplus} x)}

  where \eqn{x} is a unit vector in \eqn{R^q}, and \code{mu} (\eqn{q}-dimensional vector) and
  \code{Aplus} (\eqn{q} by \eqn{q} symmetric) are parameters.

  The matrix Bingham distribution on \eqn{q} by \eqn{r} matrices \eqn{X} whose columns are
  orthonormal, is given by the density proportional to 

  \deqn{\exp(trace(X^T \code{Aplus} X)).}

  The Bingham distribution on the unit sphere in \eqn{R^q} can be simulated
  using (a) \code{rBingham}, (b) \code{rFisherBingham} with \code{mu=0}, and (c)
  \code{rBingham.Grassmann} with \code{r=1}.  Choice (a) will be fastest.

  The Fisher distribution can be simulated using \code{rFisherBingham} with
  \code{Aplus=0}.

  The matrix Fisher distribution on SO(3) has density proportional to 

  \deqn{\exp(trace(\code{Fmat}^T X))}

  where X is a 3 by 3 rotation matrix, and \code{Fmat} is a 3 by 3 parameter
  matrix.

  The bivariate von Mises sine model on the torus has density
  proportional to

  \deqn{\exp(\code{k1} \cos(\theta) + \code{k2} cos(\phi) + \code{alpha} sin(\theta) \sin(\phi))}

  for two angles \eqn{\theta}, \eqn{\phi}.  The Bessel density is obtained from the
  bivariate von Mises sine model as the marginal density of \eqn{\theta}.}

\value{For \code{rBingham} and \code{rFisherBingham}, the output is an \code{nsim} by \code{q}
  matrix.  Each row is a simulated unit vector.

  For \code{rBingham.Grassmann}, the output is an \code{nsim} by \code{q} by \code{r} array.  For
  each value of the first index, the result is a simulated \code{q} by \code{r} matrix with
  orthonormal columns.

  For \code{rFisher.SO3},  the output is an \code{nsim} by \code{r} by \code{r} array.  For
  each value of the first index, the result is a simulated  \code{r} by \code{r} matrix
  rotation matrix.

  For \code{rvMsin.torus}, the output is an \code{nsim} by 4 matrix, with each row
  containing the simulated value of \eqn{(\cos(\theta), \sin(\theta), \cos(\phi),
  \sin(\phi))}.

  For \code{rBessel}, the output is an \code{nsim} by 2 matrix containing the marginal
  simulated values of \eqn{(\cos(\theta), \sin(\theta))}.

  In all cases, the output has an attribute \code{summary}, which is a vector
  of length 6 summarizing some details about the number of simulations
  needed in the acceptance/rejection algorithm.  The key element of this
  vector is called \code{efficiency}, a number between 0 and 1, where 1 means
  that all the simulated values from the envelope distribution have been
  accepted. The \code{success} entry is \code{0} if the simulations were not completed, often due to the number of iterations (in entry \code{mloops}) reaching \code{mtop}.}

\references{
John T. Kent, Asaad M. Ganeiber & Kanti V. Mardia (2018) A New Unified Approach for the Simulation of a Wide Class of Directional Distributions, Journal of Computational and Graphical Statistics, 27:2, 291-301, DOI: \doi{10.1080/10618600.2017.1390468} 
}
\author{J T Kent
}
\seealso{
  The function \code{\link[CircStats]{rvm}} in the CRAN library \code{CircStats} provides a more
  efficient method to simulate from the von Mises distribution (i.e. the
  Fisher distribution with q=2) than \code{rFisherBingham}.
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
