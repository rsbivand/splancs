% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{plt}
\alias{plt}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{bins nearest neighbour distances}
\description{
 bins nearest neighbour distances
}
\usage{
plt(data, value)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{nearest neighbour distances}
  \item{value}{breaks for binning distances}
}

\value{
  binned values
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

\seealso{\code{\link{Fhat}}, \code{\link{Ghat}}}

