% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{mpoint}
\alias{mpoint}
\title{
Overlay a number of point patterns
}
\description{
Overlay a number of point patterns.
}
\usage{
mpoint(\dots,cpch,add=FALSE,type="p")
}
\arguments{
\item{...}{
At least one argument consisting of a points data set must be specified.
}
\item{cpch}{
A vector of characters for plotting symbols
}
\item{add}{
if \code{add} is \code{TRUE} then overlay on an existing plot
}
\item{type}{
plot data as points if \code{type="p"}, lines if \code{type="l"}


}}
\details{
\code{mpoint} enables several point or polygon datasets to be overlayed. The plot
region is calculated so that all the specified datasets fit in the region.
The parameter \code{cpch} specifies the characters to use for each set of points. The
default \code{cpch} consists of the numbers 1 to 9 followed by the uppercase 
letters A to Z. If \code{cpch} is shorter than the number of point sets to
plot, then it is repeated.
}
\seealso{
\code{\link{pointmap}}
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

