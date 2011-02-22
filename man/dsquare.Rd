% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{dsquare}
\alias{dsquare}
\title{
Distance-squared from a number of points to a number of sources
}
\description{
Computes the distance-squared from a number of points to a number of sources.
}
\usage{
dsquare(pts, srcs, namepref="d")
}
\arguments{
\item{pts}{
A number of points representing the locations of cases and controls.
}
\item{srcs}{
A number of points representing source locations
}
\item{namepref}{
A prefix given to the name of the results.
}}
\value{
A data frame with the same number of  columns as \code{srcs}. The column names will 
be the value of \code{namepref} prefixing the numbers from 1 to the number
of sources.
}
\seealso{
\code{\link{tribble}}, \code{\link{triblik}}


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

