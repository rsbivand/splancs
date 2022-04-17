% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{kerview}
\alias{kerview}
\title{
A linked-window system for browsing space-time data
}
\description{
A linked-window system for browsing space-time data.
}
\usage{
kerview(pts, times, k3, map=TRUE, addimg=TRUE, ncol=1)
}
\arguments{
\item{pts}{
A matrix of event x,y coordinates.
}
\item{times}{
A vector of event times.
}
\item{k3}{
An object returned from \code{kernel3d}, the space-time kernel smoothing
function
}
\item{map}{
If false, don't plot the map display.
}
\item{addimg}{
If true, overwrite successive images in the image display, else
make a fresh image plot each time.
}
\item{ncol}{
Number of columns and rows for multiple images and maps.
}}
\details{
This function displays three linked views of the data. In the current
graphics window a temporal slice from the kernel smoothing is displayed. 
Another graphics device is started to display a map of the data that
contributed to that time-slice. A third graphics device shows a 
histogram of the times of the events. Clicking with the mouse
in this window with button 1 sets the time for the other
displays to the time on the x-axis of the histogram at the clicked point.


In this way the 3-dimensional kernel smoothed function can be browsed, 
and the corresponding map of the data compared.


}
\seealso{
\code{\link{kernel3d}}


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

