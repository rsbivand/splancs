% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{khvc}
\alias{khvc}
\title{
Covariance matrix for the difference between two K-functions
}
\description{
Calculate the covariance matrix for the difference between two K-functions.
Also return the contribution to the variance for each of the two point
patterns,
}
\usage{
khvc(pts1, pts2, poly, s)
}
\arguments{
\item{pts1}{
An object containing the case locations.
}
\item{pts2}{
An object containing the control locations.
}
\item{poly}{
A polygon enclosing the locations in \code{pts1} and \code{pts2}
}
\item{s}{
A vector of distances at which the calculation is to be made.
}}
\value{
A list with four components:
\item{varmat}{The upper triangle of the covariance matrix.}
\item{k11}{The variance of Khat for the cases}
\item{k22}{The variance of Khat for the controls}
\item{k12}{The covariance of Khat for the cases and Khat for controls.}
}
\section{Note}{
Note that the diagonal of the covariance matrix is
 \code{$k11 - 2 * $k12 + $k22
}}
\seealso{
\code{\link{khat}}, \code{\link{khvmat}}, \code{\link{secal}}


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

