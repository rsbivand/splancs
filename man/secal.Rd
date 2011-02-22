% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{secal}
\alias{secal}
\title{
Standard errors for the difference between two K-functions
}
\description{
Calculate standard errors for the difference between two K-functions under
random labelling of the corresponding two sets of points.
}
\usage{
secal(pts1,pts2,poly,s)
}
\arguments{
\item{pts1,pts2}{
Two point data sets
}
\item{poly}{
Polygon enclosing the points in \code{pts1} and \code{pts2}
}
\item{s}{
A vector of distances at which to calculate the standard error.
}}
\value{
A vector like \code{s} containing the value of the standard error at each of the
distances in \code{s}


}
\details{
To compare two point patterns, one can calculate the difference between
their K-functions. The function \code{secal} gives the pointwise standard 
errors for the estimated differences, under the random labelling hypothesis.
}
\references{
Diggle P.J. and Chetwynd A.G. (1991) Second-order analysis of spatial 
clustering
\emph{Biometrics 47} 1155--63; 
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{http://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

\seealso{
\code{\link{khat}}
}
% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

