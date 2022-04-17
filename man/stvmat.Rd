% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{stvmat}
\alias{stvmat}
\title{
Variance matrix for space-time clustering
}
\description{
Compute the variance matrix for space-time clustering
}
\usage{
stvmat(pts, times, poly, tlim, s, tm)
}
\arguments{
\item{pts}{
A set of points.
}
\item{times}{
A vector of times, the same length as the number of points in \code{pts}
}
\item{poly}{
A polygon that encloses the points
}
\item{tlim}{
A vector of length 2 specifying the upper and lower temporal domain.
}
\item{s}{
A vector of spatial distances for the analysis
}
\item{tm}{
A vector of times for the analysis
}}
\value{
A four-dimensional matrix is returned. The covariance between space-time
\code{t1,s1} and \code{t2,s2} is given by the corresponding element \code{[t1,s1,t2,s2]}
For full details, see Diggle, Chetwynd, Haggkvist and Morris (1995)
}
\seealso{
\code{\link{stkhat}}, \code{\link{stsecal}}, \code{\link{stmctest}}, \code{\link{stdiagn}}
}

\references{
Diggle, P., Chetwynd, A., Haggkvist, R. and Morris, S. 1995 Second-order analysis of space-time clustering. Statistical Methods in Medical Research, 4, 124-136;
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{https://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

