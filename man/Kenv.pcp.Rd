% Copyright Giovanni Petris and Roger S. Bivand 2001
%

\name{Kenv.pcp}
\alias{Kenv.pcp}
\title{Calculate simulation envelope for a Poisson Cluster Process}
\description{
  This function computes the envelope of Khat from simulations of a Poisson Cluster Process for a given polygon
}
\usage{
Kenv.pcp(rho, m, s2, region.poly, larger.region=NULL, nsim, r, vectorise.loop=TRUE)
}
\arguments{
  \item{rho}{intensity of the parent process}
  \item{m}{average number of offsprings per parent}
  \item{s2}{variance of location of offsprings relative to their parent}
  \item{region.poly}{a polygon defining the region in which the process is to be generated}
  \item{larger.region}{a rectangle containing the region of interest given in the form (xl,xu,yl,yu), defaults to \code{sbox()} around region.poly}
  \item{nsim}{number of simulations required}
  \item{r}{vector of distances at which the K function has to be estimated}
  \item{vectorise.loop}{if TRUE, use new vectorised code, if FALSE, use loop as before}
}
\value{
  \item{ave}{mean of simulations}
  \item{upper}{upper bound of envelope}
  \item{lower}{lower bound of envelope}
}
\references{Diggle, P. J. (1983) \emph{Statistical analysis of spatial point
patterns}, London: Academic Press, pp. 55-57 and 78-81; Bailey, T. C. and
Gatrell, A. C. (1995) \emph{Interactive spatial data analysis}, Harlow:
Longman, pp. 106-109.}
\author{Giovanni Petris <GPetris@uark.edu>,  Roger.Bivand@nhh.no}

\seealso{\code{\link{pcp}}, \code{\link{pcp.sim}}, \code{\link{khat}}}

\examples{
data(cardiff)
polymap(cardiff$poly)
pointmap(as.points(cardiff), add=TRUE)
title("Locations of homes of 168 juvenile offenders")
pcp.fit <- pcp(as.points(cardiff), cardiff$poly, h0=30, n.int=30)
pcp.fit
m <- npts(as.points(cardiff))/(areapl(cardiff$poly)*pcp.fit$par[2])
r <- seq(2,30,by=2)
K.env <- Kenv.pcp(pcp.fit$par[2], m, pcp.fit$par[1], cardiff$poly,
           nsim=20, r=r)
L.env <- lapply(K.env, FUN=function(x) sqrt(x/pi)-r)
limits <- range(unlist(L.env))
plot(r, sqrt(khat(as.points(cardiff),cardiff$poly,r)/pi)-r, ylim=limits,
     main="L function with simulation envelopes and average", type="l",
     xlab="distance", ylab="")
lines(r, L.env$lower, lty=5)
lines(r, L.env$upper, lty=5)
lines(r, L.env$ave, lty=6)
abline(h=0)
}
\keyword{spatial}
