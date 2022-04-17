% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{Kenv.csr}
\alias{Kenv.csr}
\title{
Envelope of Khat from simulations of complete spatial randomness
}
\description{
Compute envelope of Khat from simulations of complete spatial randomness.
}
\usage{
Kenv.csr(nptg,poly,nsim,s,quiet=FALSE)
}
\arguments{
\item{nptg}{
Number of points to generate in each simulation.
}
\item{poly}{
Polygon in which to generate the points.
}
\item{nsim}{
Number of simulations to do.
}
\item{s}{
Vector of distances at which to calculate the envelope.
}
\item{quiet}{
If FALSE, print a message after every simulation for progress monitoring. If
TRUE, print no messages.
}}
\value{
A list with two components, called \code{$upper} and \code{$lower}. Each
component is a vector like \code{s}. The two components contain the upper 
and lower bound of the Khat envelope. 


}
\examples{
data(cardiff)
UL.khat <- Kenv.csr(length(cardiff$x), cardiff$poly, nsim=29, seq(2,30,2))
plot(seq(2,30,2), sqrt(khat(as.points(cardiff), cardiff$poly, 
seq(2,30,2))/pi)-seq(2,30,2), type="l", xlab="Splancs - polygon boundary", 
ylab="Estimated L", ylim=c(-1,1.5))
lines(seq(2,30,2), sqrt(UL.khat$upper/pi)-seq(2,30,2), lty=2)
lines(seq(2,30,2), sqrt(UL.khat$lower/pi)-seq(2,30,2), lty=2)
}
\seealso{
\code{\link{csr}}, \code{\link{khat}}
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

