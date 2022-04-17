% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2001 by Roger S. Bivand

\name{inout}
\alias{inout}
\title{
Test points for inclusion in a polygon
}
\description{
Test points for inclusion in a polygon.
}
\usage{
inout(pts,poly,bound=NULL,quiet=TRUE)
}
\arguments{
\item{pts}{
A points data set
}
\item{poly}{
A polygon data set
}
\item{bound}{If points fall exactly on polygon boundaries, the default NULL gives arbitrary assignments. If TRUE, then all points "on" boundaries are set as within the polygon, if FALSE, outside.}
\item{quiet}{Do not report which points are on boundary for non-NULL bound}
}
\value{
A vector of logical values. \code{TRUE} means the point was inside the
polygon, \code{FALSE} means the point was outside. Note that "inside" is an arbitrary concept for points "on" the polygon boundary.



}
\seealso{
\code{\link{inpip}},\code{\link{pip}}
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
data(uganda)
suganda <- sbox(uganda$poly)
ruganda <- csr(suganda, 1000)
polymap(suganda)
polymap(uganda$poly, add=TRUE)
def <- inout(ruganda, uganda$poly, bound=NULL)
pointmap(as.points(ruganda[def,1], ruganda[def,2]), add=TRUE, col="black")
pointmap(as.points(ruganda[!def,1], ruganda[!def,2]), add=TRUE, col="red")
tru <- inout(ruganda, uganda$poly, bound=TRUE, quiet=FALSE)
which(tru & !def)
ds1 <- as.points(expand.grid(x=seq(-1.5,1.5,0.5), y=seq(-1.5,1.5,0.5)))
ds1.poly <- ds1[chull(ds1),]
ds2 <- as.points(rnorm(300),rnorm(300))
plot(ds2, type="n", asp=1)
polymap(ds1.poly, add=TRUE, border="lightblue", col="lightblue", lwd=1)
points(ds2[inout(ds2,ds1.poly),],  col="green",  pch=20)
points(ds2[!inout(ds2,ds1.poly),], col="orange", pch=20)
points(ds1[inout(ds1,ds1.poly),],  col="black",  pch=20)
points(ds1[!inout(ds1,ds1.poly),], col="red",    pch=20)
plot(ds2, type="n", asp=1)
polymap(ds1.poly, add=TRUE, border="lightblue", col="lightblue", lwd=1)
points(ds2[inout(ds2,ds1.poly,bound=TRUE),],  col="green",  pch=20)
points(ds2[!inout(ds2,ds1.poly,bound=TRUE),], col="orange", pch=20)
points(ds1[inout(ds1,ds1.poly,bound=TRUE),],  col="black",  pch=20)
points(ds1[!inout(ds1,ds1.poly,bound=TRUE),], col="red",    pch=20)
plot(ds2, type="n", asp=1)
polymap(ds1.poly, add=TRUE, border="lightblue", col="lightblue", lwd=1)
points(ds2[inout(ds2,ds1.poly,bound=FALSE),],  col="green",  pch=20)
points(ds2[!inout(ds2,ds1.poly,bound=FALSE),], col="orange", pch=20)
points(ds1[inout(ds1,ds1.poly,bound=FALSE),],  col="black",  pch=20)
points(ds1[!inout(ds1,ds1.poly,bound=FALSE),], col="red",    pch=20)
}

\keyword{spatial}

