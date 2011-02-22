% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{getpoly}
\alias{getpoly}
\title{
Draw a polygon on the current graphics device
}
\description{
Draw a polygon on the current graphics device
}
\usage{
getpoly(quiet=FALSE)
}
\arguments{
\item{quiet}{
if TRUE, don't prompt for input of a polygon.
}}
\value{
A polygon data set consisting of the points entered. The current coordinate 
system is used.


}
\details{
The system prompts the user to enter points on the current graphics
device using the mouse or other pointing device. The points are joined
on the screen with the current line symbol.
A polygon of the points entered is drawn on the current graphics device.
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

