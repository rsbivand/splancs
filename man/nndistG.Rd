% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{nndistG}
\alias{nndistG}
\title{
Nearest neighbour distances as used by Ghat()
}
\description{
Calculate nearest neighbour distances as used by \code{Ghat()}.
}
\usage{
nndistG(pts)
}
\arguments{
\item{pts}{
A points data set
}}
\value{
Returns a list with components \code{$dists} and \code{$neighs}.
\code{$dists[i]} is the distance
to the nearest neighbour of point \code{i} in \code{pts,} and \code{$neighs[i]}
is the index
of the neighbour of point \code{i}.


}
\examples{
data(uganda)
boxplot(nndistG(as.points(uganda))$dists)
plot(ecdf(nndistG(as.points(uganda))$dists))
}
\seealso{
\code{\link{nndistF}}, \code{\link{Fhat}}, \code{\link{Ghat}}, \code{\link{Fzero}}
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

