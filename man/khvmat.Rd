% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{khvmat}
\alias{khvmat}
\title{
Covariance matrix for the difference between two K-functions
}
\description{
Calculate the covariance matrix for the difference between two K-functions 
under random labelling of the corresponding two sets of points.
}
\usage{
khvmat(pts1, pts2, poly, s)
}
\arguments{
\item{pts1}{
An object containing the case locations.
}
\item{pts2}{
An object containing the control locations.
}
\item{poly}{
Polygon enclosing the points in pts1 and pts2.
}
\item{s}{
A vector of distances at which the calculation is to be made.
}}
\value{
A matrix containing the covariances, with the variances on the diagonal.
}

\seealso{
\code{\link{secal}}
}

\references{Diggle P.J and Chetwynd A.C (1991) Second order analysis of
spatial clustering Biometrics 47 1155-63;
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{http://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

% Converted by Sd2Rd version 0.3-3.

\keyword{spatial}

