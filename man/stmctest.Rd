% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2009 by Roger S. Bivand

\name{stmctest}
\alias{stmctest}
\title{
Monte-Carlo test of space-time clustering
}
\description{
Perform a Monte-Carlo test of space-time clustering.
}
\usage{
stmctest(pts, times, poly, tlimits, s, tt, nsim, quiet=FALSE, returnSims=FALSE)
}
\arguments{
\item{pts}{
A set of points as used by Splancs.
}
\item{times}{
A vector of times, the same length as the number of points in \code{pts}.
}
\item{poly}{
A polygon enclosing the points.
}
\item{tlimits}{
A vector of length 2, specifying the upper and lower temporal domain.
}
\item{s}{
A vector of spatial distances for the analysis.
}
\item{tt}{
A vector of times for the analysis.
}
\item{nsim}{
The number of simulations to do.
}
\item{quiet}{
If \code{quiet=TRUE} then no output is produced, otherwise the function
prints the number of simulations completed so far, and also how the
test statistic for the data ranks with the simulations.
}
\item{returnSims}{default FALSE, if TRUE, return the \code{stkhat} output for the observed data and each simulation as attributes \code{obs} and \code{sims}}
}
\value{
A list with components:

\item{t0}{The observed value of the statistic}
\item{t}{A single column matrix with \code{nsim} values each of which is a simulated value of the statistic}
%A vector of length \code{nsim+1}. The first element is the test statistic
%for the data, and the remaining elements are those for the simulations.
}
\details{
The function  uses a sum of residuals as a test statistic, 
randomly permutes the times of the set of points and recomputes the
test statistic for a number of simulations.
See Diggle, Chetwynd, Haggkvist and Morris (1995) for details.
}
\seealso{
\code{\link{stkhat}}, \code{\link{stsecal}}, \code{\link{stvmat}}, \code{\link{stdiagn}}
}

\references{
Diggle, P., Chetwynd, A., Haggkvist, R. and Morris, S. 1995 Second-order analysis of space-time clustering. Statistical Methods in Medical Research, 4, 124-136;Bailey, T. C. and Gatrell, A. C. 1995, Interactive spatial data analysis. Longman, Harlow, pp. 122-125;
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{http://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

\note{The example of using returned simulated values is included only to show how the values might be used, not to indicate that this constitutes a way of examining which observed values of the space-time measure are exceptional.}

\examples{
example(stkhat)
bur1mc <- stmctest(burpts, burkitt$t, burbdy, c(400, 5800),
  seq(1,40,2), seq(100, 1500, 100), nsim=49, quiet=TRUE, returnSims=TRUE)
plot(density(bur1mc$t), xlim=range(c(bur1mc$t0, bur1mc$t)))
abline(v=bur1mc$t0)
r0 <- attr(bur1mc, "obs")$kst-outer(attr(bur1mc, "obs")$ks, attr(bur1mc, "obs")$kt)
rsimlist <- lapply(attr(bur1mc, "sims"), function(x) x$kst - outer(x$ks, x$kt))
rarray <- array(do.call("cbind", rsimlist), dim=c(20, 15, 49))
rmin <- apply(rarray, c(1,2), min)
rmax <- apply(rarray, c(1,2), max)
r0 < rmin
r0 > rmax
}
% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

