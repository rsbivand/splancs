% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{Shift}
\alias{Shift}
\title{
Shift a point data set
}
\description{
Shift a point data set (function name changed from shift to Shift to avoid collision with spatstat)
}
\usage{
Shift(pts,xsh=0.0,ysh=0.0)
}
\arguments{
\item{pts}{
The point data set to shift
}
\item{xsh}{
Amount to shift along the x-axis
}
\item{ysh}{
Amount to shift along the y-axis
}}
\value{
A point data set like \code{pts}, but with \code{xsh} added to its
x-coordinates, and \code{ysh} added to its y-coordinates.


}
\seealso{
\code{\link{tor.shift}}
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

