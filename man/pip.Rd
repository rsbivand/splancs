% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{pip}
\alias{pip}
\title{
Points inside or outside a polygon
}
\description{
Return points inside or outside a polygon.
}
\usage{
pip(pts,poly,out=FALSE,bound=NULL,quiet=TRUE)
}
\arguments{
\item{pts}{
A points data set
}
\item{poly}{
A polygon data set
}
\item{out}{
If \code{out=TRUE}, return the points outside the polygon, else the points inside.
}
\item{bound}{If points fall exactly on polygon boundaries, the default NULL gives arbitrary assignments. If TRUE, then all points "on" boundaries are set as within the polygon, if FALSE, outside.}
\item{quiet}{Do not report which points are on boundary for non-NULL bound}
}
\value{
\code{pip} returns the points of \code{pts} that lie inside (or outside with
\code{out=TRUE})
the polygon \code{poly}. Compare this with \code{inpip}, which returns
the indices of
the points in the polygon, and \code{inout} which returns a logical vector
indicating whether points are inside or outside.


}
\details{
\code{pip} calls \code{inout}, then selects the appropriate sub-set of points.
}
\seealso{
\code{\link{inpip}}, \code{\link{inout}}
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

