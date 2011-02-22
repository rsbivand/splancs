% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2001 by Roger S. Bivand

\name{inpip}
\alias{inpip}
\title{
Select points inside a polygon
}
\description{
Select points inside a polygon
}
\usage{
inpip(pts,poly,bound=NULL,quiet=TRUE)
}
\arguments{
\item{pts}{
A points data set
}
\item{poly}{
A polygon data set
}
\item{bound}{If points fall exactly on polygon boundaries, the default NULL gives arbitrary assignments. If TRUE, then all points "on" boundaries are set as within the polygon, if FALSE, outside.}
\item{quiet}{Do not report which points are on boundary for non-NULL bound}
}
\value{
\code{inpip} returns a vector of indices of the points in \code{pts} that are located
in the polygon. Note that "in" is an arbitrary concept for points "on" the polygon boundary.


}
\seealso{
\code{\link{pip}}
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

