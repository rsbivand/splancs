% Copyright Giovanni Petris and Roger S. Bivand 2001

\name{pcp}
\alias{pcp}
\title{Fit a Poisson cluster process}
\description{
  The function fits a Poisson cluster process to point data for a given
enclosing polygon and fit parameters
}
\usage{
pcp(point.data, poly.data, h0=NULL, expo=0.25, n.int=20)
}
\arguments{
  \item{point.data}{a points object}
  \item{poly.data}{a polygon enclosing the study region}
  \item{h0}{upper bound of integration in the criterion function}
  \item{expo}{exponent in the criterion function}
  \item{n.int}{number of intervals used to approximate the integral in the criterion function with a sum}
}
\value{
  The function returns an object as returned by \code{optim}, including:
  \item{par}{The best set of parameters s2 and rho found}
  \item{value}{The value of the fit corresponding to `par'}
  \item{convergence}{`0' indicates successful convergence}
}
\references{Diggle, P. J. (1983) \emph{Statistical analysis of spatial point
patterns}, London: Academic Press, pp. 55-57 and 78-81; Bailey, T. C. and
Gatrell, A. C. (1995) \emph{Interactive spatial data analysis}, Harlow:
Longman, pp. 106-109.}
\author{Giovanni Petris <GPetris@uark.edu>,  Roger.Bivand@nhh.no}

\seealso{\code{\link{optim}}, \code{\link{pcp.sim}}, \code{\link{Kenv.pcp}}, \code{\link{khat}}}

\examples{
data(cardiff)
polymap(cardiff$poly)
pointmap(as.points(cardiff), add=TRUE)
title("Locations of homes of 168 juvenile offenders")
pcp.fit <- pcp(as.points(cardiff), cardiff$poly, h0=30, n.int=30)
pcp.fit
}
\keyword{spatial}
