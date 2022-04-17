% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{areapl}
\alias{areapl}
\title{
Calculate area of polygon
}
\description{
Calculate area of polygon. If the polygon is self-intersecting, the area will not be correct.
}
\usage{
areapl(poly)
}
\arguments{
\item{poly}{
a polygon data set
}}
\value{
The area of the polygon is returned


}

\references{
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{https://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

\examples{
x <- c(1,0,0,1,1,1,1,3,3,1)
y <- c(0,0,1,1,0,0,-1,-1,0,0)
m <- cbind(x, y)
plot(m, type="b")
areapl(m)
areapl(m[1:5,])
areapl(m[6:10,])
}

% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

