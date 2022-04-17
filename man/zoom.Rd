% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{zoom}
\alias{zoom}
\title{
Interactively specify a region of a plot for expansion
}
\description{
Interactively specify a region of a plot for expansion
}
\usage{
zoom(quiet=FALSE,out=FALSE,...)
}
\arguments{
\item{quiet}{
If false, prompt the user to enter two coordinates. If true, say nothing.
}
\item{out}{
If true, expand the limits of the current plot by a factor of three, centred
on the current plot.
}
\item{...}{
Other arguments are passed through to pointmap.
}}
\value{
None


}
\details{
A prompt is optionally displayed, and the user selects two points
forming the diagonal of a rectangle. A new, empty plot is created that has its
axis limits set to the bounding square of the selected rectangle. 
If \code{out=TRUE}, no prompt is displayed, and a new blank plot is created with
its limits in x and y set to span an area three times the height and width
centred on the current centre.
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

