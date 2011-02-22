% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{spoints}
\alias{spoints}
\title{
Point Objects
}
\description{
Creates and tests for data in spatial point format.
}
\usage{
spoints(data,npoints)
}
\arguments{
\item{data}{
vector containing the data values for the points
in order (x1,y1),(x2,y2),\dots{}


}
\item{npoints}{
number of points to generate, if missing, set to length(data)/2.
}}
\value{
\code{spoints} returns an object suitable for use as a point data
object. If \code{npoints} is given, the vector data is either truncated
or repeated until sufficient data values are generated.
The returned object is a two-column matrix, where the first column stores the 
x-coordinate, and the second column stores the y-coordinate.


}
\seealso{
\code{\link{npts}}
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

