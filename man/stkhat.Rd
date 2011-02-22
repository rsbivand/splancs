% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{stkhat}
\alias{stkhat}
\title{
Space-time K-functions
}
\description{
Compute the space-time K-functions
}
\usage{
stkhat(pts, times, poly, tlimits, s, tm)
}
\arguments{
\item{pts}{
A set of points as defined in Splancs
}
\item{times}{
A vector of times, the same length as the number of points in \code{pts}
}
\item{poly}{
A polygon enclosing the points
}
\item{tlimits}{
A vector of length 2 specifying the upper and lower temporal domain.
}
\item{s}{
A vector of spatial distances for the analysis.
}
\item{tm}{
A vector of times for the analysis
}}
\value{
A list with the following components is returned:
\item{s, t}{The spatial and temporal scales}
\item{ks}{The spatial K-function}
\item{kt}{The temporal K-function}
\item{kst}{The space-time K-function}
For details see Diggle, Chetwynd, Haggkvist, and Morris (1995)
}
\seealso{
\code{\link{stsecal}}, \code{\link{stvmat}}, \code{\link{stmctest}}, \code{\link{stdiagn}}
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
\examples{
data(burkitt)
bur1 <- stkhat(burpts, burkitt$t, burbdy, c(400, 5800),
  seq(1,40,2), seq(100, 1500, 100))
oldpar <- par(mfrow=c(2,1))
plot(bur1$s, bur1$ks, type="l", xlab="distance", ylab="Estimated K",
  main="spatial K function")
plot(bur1$t, bur1$kt, type="l", xlab="time", ylab="Estimated K",
  main="temporal K function")
par(oldpar)
}
% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

