% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{sbox}
\alias{sbox}
\title{
Generate a box surrounding a point object
}
\description{
Generate a box surrounding a point object
}
\usage{
sbox(pts, xfrac = .1, yfrac = .1)
}
\arguments{
\item{pts}{
A points data set
}
\item{xfrac}{
The fraction of the width of the point pattern by which the box will
surround the point pattern to the left and right.
}
\item{yfrac}{
The fraction of the height of the point pattern by which the box will 
surround the point pattern to the top and bottom. 
}}
\value{
A points data set of four points giving the coordinates of the surrounding box


}
\seealso{
\code{\link{bboxx}}
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

