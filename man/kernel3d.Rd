% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{kernel3d}
\alias{kernel3d}
\title{
Space-time kernel
}
\description{
Compute the space-time kernel
}
\usage{
kernel3d(pts, times, xgr, ygr, zgr, hxy, hz)
}
\arguments{
\item{pts}{
A matrix of event coodinates x,y.
}
\item{times}{
A vector of event times, t.
}
\item{xgr}{
The values of x at which to compute the kernel function.
}
\item{ygr}{
The values of y at which to compute the kernel function.
}
\item{zgr}{
The values of time at which to compute the kernel function.
}
\item{hxy}{
The quartic kernel width in the x and y direction.
}
\item{hz}{
The quartic kernel width in the temporal direction.
}}
\value{
A list is returned. Most of the components are just copies of the 
input parameters, except for the \code{$v} parameter.
This is a three dimensional array containing the kernel-smoothed
values. Its dimension is
\code{[length(xgr),length(ygr),length(tgr)]}. 
}
\seealso{
\code{\link{kerview}}
}

\references{
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{https://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

\examples{
data(burkitt)
b3d <- kernel3d(burpts, burkitt$t, seq(250,350,10), seq(250, 400, 10),
  seq(365,5800,365), 30, 200)
brks <- quantile(b3d$v, seq(0,1,0.05))
cols <- heat.colors(length(brks)-1)
oldpar <- par(mfrow=c(3,5))
for (i in 1:15) image(seq(250,350,10), seq(250, 400, 10), b3d$v[,,i],
  asp=1, xlab="", ylab="", main=1960+i, breaks=brks, col=cols)
par(oldpar)
}

% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

