% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{Ghat}
\alias{Ghat}
\title{
G nearest neighbour distribution function
}
\description{
Calculates an estimate of the G nearest neighbour distribution function.
}
\usage{
Ghat(pts,s)
}
\arguments{
\item{pts}{
A points data set
}
\item{s}{
A vector of distances at which to evaluate the G function
}}
\value{
A vector of the same length as \code{s}, containing the estimate of G at the 
distances in \code{s}.


}
\details{
The function \code{Ghat(pts,s)} is defined as the proportion of members of a
point set for which the distance to the nearest other member of the set
is less than or equal to \code{s}.
}
\examples{
data(uganda)
plot(seq(20, 500, 20), Ghat(as.points(uganda), seq(20, 500, 20)), 
type="l", xlab="distance", ylab="Estimated G")
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
\seealso{\code{\link{Fhat}}, \code{\link{nndistG}}}
