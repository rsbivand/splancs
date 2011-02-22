% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{Fhat}
\alias{Fhat}
\title{
F nearest neighbour distribution function
}
\description{
Calculates an estimate of the F nearest neighbour distribution function
}
\usage{
Fhat(pts1,pts2,s)
}
\arguments{
\item{pts1}{
A points data set
}
\item{pts2}{
A points data set
}
\item{s}{
A vector of distances at which to evaluate Fhat
}}
\value{
A vector of the same length as \code{s}, containing the value of Fhat at the 
distances in \code{s}.


}
\details{
The function \code{Fhat(pts1,pts2,s)} is defined as the proportion of members of a
point set \code{pts2} for which the distance to the nearest member of
another points set \code{pts1} is less than or equal to \code{s}.
}
\examples{
data(uganda)
plot(seq(20, 500, 20), Fhat(as.points(uganda), 
as.points(csr(uganda$poly, length(uganda$x))), seq(20, 500, 20)), 
type="l", xlab="distance", ylab="Estimated F")
plot(Ghat(as.points(uganda), seq(20, 500, 20)), Fhat(as.points(uganda), 
as.points(csr(uganda$poly, length(uganda$x))), seq(20, 500, 20)), 
type="l", xlab="Estimated G", ylab="Estimated F")
lines(c(0,1),c(0,1),lty=2)
}
\seealso{
\code{\link{nndistF}},\code{\link{Fzero}}
}

\references{
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{http://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

