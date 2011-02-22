% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{bboxx}
\alias{bboxx}
\title{
Generate a non-closed bounding polygon
}
\description{
Generate a non-closed bounding polygon from the bounding box of an object
}
\usage{
bboxx(obj)
}
\arguments{
\item{obj}{
A matrix with two rows and two columns reporting the bounding box of an object
}}

\details{
The object used by \code{bboxx} may easily be created by using the \pkg{sp} \code{bbox} method on an object of interest, such as a points data set.
}

\value{
A points data set of four points giving the non-closed coordinates of the bounding box
}
\seealso{
\code{\link{sbox}}, \code{\link[sp]{bbox}}
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

