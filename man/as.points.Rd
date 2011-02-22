% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{as.points}
\alias{as.points}
\title{
Creates data in spatial point format
}
\description{
Creates data in spatial point format.
}
\usage{
as.points(...)
}
\arguments{
\item{...}{
any object(s), such as x and y vectors of the same length, or a list or data
frame containing x and y vectors. Valid options for \code{...} are:
a points object ; returns it unaltered;
a list with x and y elements of the same length --- returns a
points object with the x and y elements as the coordinates of the points;
two vectors of equal length ; returns a points object with the 
first vector as the x coordinates, the second vector as the y-coordinates.
}
}

\value{
\code{as.points} tries to return the argument(s) as a points object. }

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

