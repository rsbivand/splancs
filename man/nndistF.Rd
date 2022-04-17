% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{nndistF}
\alias{nndistF}
\title{
Nearest neighbour distances as used by Fhat()
}
\description{
Calculate nearest neighbour distances as used by  \code{Fhat()}
}
\usage{
nndistF(pts1,pts2)
}
\arguments{
\item{pts1}{
A points data set
}
\item{pts2}{
A points data set
}}
\value{
The set of distances from each of the points in \code{pts2} to the nearest
point in \code{pts1} is returned as a vector.


}
\seealso{
\code{\link{nndistG}}, \code{\link{Fhat}}, \code{\link{Ghat}}, \code{\link{Fzero}}
}
\examples{
data(uganda)
boxplot(nndistF(as.points(uganda), as.points(csr(uganda$poly, length(uganda$x)))))
plot(ecdf(nndistF(as.points(uganda), 
as.points(csr(uganda$poly, length(uganda$x))))),
main="Fhat ecdf Uganda volcano data")
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

