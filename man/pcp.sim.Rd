% Copyright Giovanni Petris and Roger S. Bivand 2001

\name{pcp.sim}
\alias{pcp.sim}
\title{Generate a Poisson Cluster Process}
\description{
  The function generates a Poisson cluster process for a given polygon within
a larger bounding region and given process parameters
}
\usage{
pcp.sim(rho, m, s2, region.poly, larger.region=NULL, vectorise.loop=TRUE)
}
\arguments{
  \item{rho}{intensity of the parent process}
  \item{m}{average number of offsprings per parent}
  \item{s2}{variance of location of offsprings relative to their parent}
  \item{region.poly}{a polygon defining the region in which the process is to be generated}
  \item{larger.region}{a rectangle containing the region of interest given in the form (xl,xu,yl,yu), defaults to \code{sbox()} around region.poly}
  \item{vectorise.loop}{if TRUE, use new vectorised code, if FALSE, use loop as before}
}
\details{
  The function generates the parents in the larger bounding region, generates
their children also in the larger bounding region, and then returns those
inside the given polygon.}
\value{
  A point object with the simulated pattern
}
\references{Diggle, P. J. (1983) \emph{Statistical analysis of spatial point
patterns}, London: Academic Press, pp. 55-57 and 78-81; Bailey, T. C. and
Gatrell, A. C. (1995) \emph{Interactive spatial data analysis}, Harlow:
Longman, pp. 106-109.}
\author{Giovanni Petris <GPetris@uark.edu>,  Roger.Bivand@nhh.no}

\seealso{\code{\link{pcp}}, \code{\link{Kenv.pcp}}, \code{\link{khat}}}

\examples{
data(cardiff)
polymap(cardiff$poly)
pointmap(as.points(cardiff), add=TRUE)
title("Locations of homes of 168 juvenile offenders")
pcp.fit <- pcp(as.points(cardiff), cardiff$poly, h0=30, n.int=30)
pcp.fit
m <- npts(as.points(cardiff))/(areapl(cardiff$poly)*pcp.fit$par[2])
sims <- pcp.sim(pcp.fit$par[2], m, pcp.fit$par[1], cardiff$poly)
pointmap(as.points(sims), add=TRUE, col="red")
}
\keyword{spatial}
