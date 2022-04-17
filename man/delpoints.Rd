% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{delpoints}
\alias{delpoints}
\title{
Select points to delete from a points data set
}
\description{
Select points to delete from a points data set.
}
\usage{
delpoints(pts,add=FALSE)
}
\arguments{
\item{pts}{
a points data set
}
\item{add}{
if false, plot the points using \code{pointmap}.
}}
\value{
A points object containing the undeleted points.


}
\details{
Using the mouse, the user selects points on the current graphics device.
These points are marked on the plot as they are selected. The function
returns the remaining points as a points object.
If \code{add} is false the points are plotted on the current plot device.
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

