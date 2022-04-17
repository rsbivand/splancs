% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{Kenv.tor}
\alias{Kenv.tor}
\title{
Envelope of K12hat from random toroidal shifts of two point patterns
}
\description{
Compute envelope of K12hat from random toroidal shifts of two point patterns.
}
\usage{
Kenv.tor(pts1,pts2,poly,nsim,s,quiet=FALSE)
}
\arguments{
\item{pts1}{
First point data set.
}
\item{pts2}{
Second point data set.
}
\item{poly}{
Polygon containing the points.
}
\item{nsim}{
Number of random toroidal shifts to do.
}
\item{s}{
Vector of distances at which to calculate the envelope.
}
\item{quiet}{
If FALSE, print a message after every simulation for progress monitoring. If
true, print no messages.
}}
\value{
A list with two components, called \code{$upper} and \code{$lower}. Each
component is a vector like \code{s}.
}
\details{
The second point data set is randomly shifted using \code{rtor.shift}
in the rectangle defined by \code{poly}. Then \code{k12hat} is called 
to compute K12hat for the two patterns. 
The upper and lower values of K12hat over the \code{ntor}
toroidal shifts are returned.


}
\examples{
data(okwhite)
data(okblack)
okpoly <- list(x=c(okwhite$x, okblack$x), y=c(okwhite$y, okblack$y))
plot(seq(5,80,5), sqrt(k12hat(as.points(okwhite), as.points(okblack), 
bboxx(bbox(as.points(okpoly))), seq(5,80,5))/pi) - seq(5,80,5), xlab="distance", 
ylab=expression(hat(L)[12]), ylim=c(-35,35), type="l",
main="Simulation envelopes, random toroidal shifts")
env.ok <- Kenv.tor(as.points(okwhite), as.points(okblack), 
bboxx(bbox(as.points(okpoly))), nsim=29, s=seq(5,80,5))
lines(seq(5,80,5), sqrt(env.ok$upper/pi)-seq(5,80,5), lty=2)
lines(seq(5,80,5), sqrt(env.ok$lower/pi)-seq(5,80,5), lty=2)
}
\seealso{
\code{\link{rtor.shift}},\code{\link{k12hat}}
}

\references{
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{https://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

