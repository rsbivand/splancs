% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{kernrat}
\alias{kernrat}
\title{
Ratio of two kernel smoothings
}
\description{
Return the ratio of two kernel smoothings
}
\usage{
kernrat(pts1,pts2,poly,h1,h2,nx=20,ny=20,kernel='quartic')
}
\arguments{
\item{pts1,pts2}{
Point data sets
}
\item{poly}{
A polygon data set
}
\item{h1,h2}{
The kernel width parameters, \code{h1} for \code{pts1}, and \code{h2} for \code{pts2}
}
\item{nx}{
Number of points along the x-axis of the returned grid.
}
\item{ny}{
Number of points along the y-axis of the returned grid.
}
\item{kernel}{
Type of kernel function to use. Currently only the quartic kernel is implemented.
}}
\value{
A list with the following components:
\item{x}{List of x-coordinates at which the kernel function has been evaluated.}
\item{y}{List of y-coordinates at which the kernel function has been evaluated.}
\item{z}{A matrix of dimension \code{nx} by \code{ny} containing the ratio of the kernel functions.}
\item{h}{A vector of length 2 containing \code{h1} and \code{h2}}
\item{kernel}{a character string containing the kernel name.}


}
\seealso{
\code{\link{kernel2d}}, \code{\link{mse2d}}
}
\references{
Berman M. and Diggle P.J. (1989) Estimating Weighted Integrals of the
Second-Order Intensity of Spatial Point Patterns.
\emph{J. R. Statist Soc} B51 81-92;
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{http://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

\examples{
data(okwhite)
data(okblack)
okpoly <- list(x=c(okwhite$x, okblack$x), y=c(okwhite$y, okblack$y))
kr <- kernrat(as.points(okwhite), as.points(okblack), bboxx(bbox(as.points(okpoly))),
 h1=50, h2=50)
image(kr, asp=1)
brks <- quantile(c(kr$z), seq(0,1,1/10), na.rm=TRUE)
lbrks <- formatC(brks, 3, 6, "g", " ") 
cols <- heat.colors(length(brks)-1)
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,0,1,2), 2, 2, byrow = TRUE), c(2.5,1.5), c(1,3), TRUE)
image(kr, breaks=brks, col=cols, asp=1)
plot.new()
legend(c(0,1), c(0,1), legend=paste(lbrks[-length(lbrks)], lbrks[-1], sep=":"), fill=cols, bty="n")
par(def.par)
}


% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

