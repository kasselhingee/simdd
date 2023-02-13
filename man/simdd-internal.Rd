\name{simdd-internal}
\alias{b2mf}
\alias{bfind}
\alias{dimq}
\alias{dimset}
\alias{is.so3}
\alias{is.s3}
\alias{mf2b}
\alias{so3tos3}
\alias{s3toso3}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{simdd-internal}
%%  ~~function to do ... ~~

\description{Various internal routines used by simdd.  The details are
  not
documented here yet.}
%%%\usage{
%%%b2mf(A)
%%%}
%- maybe also 'usage' for other objects documented here.
%%%\arguments{
%%%  \item{A}{
%%     ~~Describe \code{A} here~~
%%%}
%%%}
%%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%%}
%%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%%}
%%%\references{
%% ~put references to the literature/web site here ~
%%%}
\author{
J T Kent
}
%%%\note{
%%  ~~further notes~~
%%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%%}
%%%\examples{
%%%##---- Should be DIRECTLY executable !! ----
%%%##-- ==>  Define data, use random,
%%%##--	or do  help(data=index)  for the standard data sets.

%%%## The function is currently defined as
%%%function (A) 
%%%{
%%%    F = matrix(0, 3, 3)
%%%    A = A - A[1, 1] * diag(4)
%%%    F[1, 1] = -(A[3, 3] + A[4, 4] - A[2, 2])/4
%%%    F[2, 2] = -(A[2, 2] + A[4, 4] - A[3, 3])/4
%%%    F[3, 3] = -(A[2, 2] + A[3, 3] - A[4, 4])/4
%%%    F[2, 3] = (A[3, 4] - A[1, 2])/2
%%%    F[3, 2] = (A[3, 4] + A[1, 2])/2
%%%    F[1, 2] = (A[2, 3] - A[1, 4])/2
%%%    F[2, 1] = (A[2, 3] + A[1, 4])/2
%%%    F[3, 1] = (A[2, 4] - A[1, 3])/2
%%%    F[1, 3] = (A[2, 4] + A[1, 3])/2
%%%    F
%%%  }
%%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%%\keyword{ ~kwd1 }
%%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
