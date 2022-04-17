% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{gen}
\alias{gen}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{generate points in polygon}
\description{
  generates random points within a defined polygon, trying to reach \code{npoints}
points - used in \code{csr}.
}
\usage{
gen(poly, npoints)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{poly}{A polygon data set}
  \item{npoints}{The number of points to generate}
}

\value{
  returns a point object.
}
\references{
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{https://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

\keyword{spatial}

\seealso{\code{\link{csr}} }


