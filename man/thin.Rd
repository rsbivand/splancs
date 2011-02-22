% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{thin}
\alias{thin}
\title{
Randomly thin a point data set
}
\description{
Randomly thin a point data set.
}
\usage{
thin(pts,n)
}
\arguments{
\item{pts}{
a points data set.
}
\item{n}{
the number of points to return
}}
\value{
Returns a point data set consisting of \code{n} points selected randomly
from the set \code{pts}.
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

