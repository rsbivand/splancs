% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{gridpts}
\alias{gridpts}
\title{
Generate a grid of points
}
\description{
Generate a grid of points
}
\usage{
gridpts(poly,npts,xs,ys)
}
\arguments{
\item{poly}{
polygon in which to generate the points
}
\item{npts}{
approximate number of points to generate
}
\item{xs,ys}{
grid spacing in x and y 


Either \code{npts} or \code{xs} and \code{ys} must be specified. If all three
are given then \code{xs} and \code{ys} are ignored. 
}}
\value{
A points object containing a grid of points inside the polygon. If \code{npts}
is specified, then a grid spacing \code{xs} and \code{ys} will be calculated to give
approximately \code{npts} in the polygon. If \code{xs} and \code{ys} are given then these
will be used to generate a number of points in the polygon.


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

