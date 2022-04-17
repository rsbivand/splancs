% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{pdense}
\alias{pdense}
\title{
Overall density for a point pattern
}
\description{
Calculate overall density for a point pattern.
}
\usage{
pdense(pts,poly)
}
\arguments{
\item{pts}{
A points data set
}
\item{poly}{
A polygon data set
}}
\value{
The density of the points in the polygon. i.e. the number of points
per unit area.


}
\seealso{
\code{\link{Fzero}}
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

