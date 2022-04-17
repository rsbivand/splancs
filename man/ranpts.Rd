% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{ranpts}
\alias{ranpts}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{adjust number of random points in polygon}
\description{
  adjust number of random points in polygon
}
\usage{
ranpts(pts, poly, nprq)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{pts}{points object}
  \item{poly}{polygon object}
  \item{nprq}{required number of points}
}

\value{
  points object with required number of random points
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

\seealso{\code{\link{csr}}}

