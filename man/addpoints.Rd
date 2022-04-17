% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{addpoints}
\alias{addpoints}
\title{
Add points interactively to a point data set
}
\description{
Add points interactively to a point data set.
}
\usage{
addpoints(pts,plot=FALSE,quiet=FALSE)
}
\arguments{
\item{pts}{
A points data set.
}
\item{plot}{
if true, plot the \code{pts} data, using \code{pointmap}. If false,
 or if \code{pts} is missing,
don't plot the data. 
}
\item{quiet}{
if true, don't print a prompt to enter points.
}}
\value{
A points data set consisting of \code{pts} and the points entered on the current 
graphics device.


}
\details{
The points entered are displayed on the current graphics device.
}
\seealso{
\code{\link{delpoints}}
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

