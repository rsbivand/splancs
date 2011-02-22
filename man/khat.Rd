% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2002 by Roger S. Bivand

\name{khat}
\alias{khat}
\alias{print.khat}
\alias{plot.khat}
\title{
K-function
}
\description{
Calculates an estimate of the K-function
}
\usage{
khat(pts,poly,s,newstyle=FALSE,checkpoly=TRUE)
\method{print}{khat}(x, ...)
\method{plot}{khat}(x, ...)
}
\arguments{
\item{pts}{
A points data set
}
\item{poly}{
A polygon containing the points - must be a perimeter ring of points
}
\item{s}{
A vector of distances at which to calculate the K function
}
\item{newstyle}{if TRUE, the function returns a khat object}
\item{checkpoly}{if TRUE compare polygon area and polygon bounding box and convex hull areas to see whether the polygon object is malformed; may be set to FALSE if the polygon is known to be a ring of points}
\item{x}{a \code{khat} object}
\item{...}{other arguments passed to plot and print functions}
}
\value{
If \code{newstyle} is FALSE, 
a vector like \code{s} containing the value of K at the points in \code{s}.
else a \code{khat} object list with:
\item{khat}{the value of K at the points in \code{s}}
\item{counts}{integer matrix of counts of points within the vector of
 distances \code{s} for each point}
\item{khats}{matrix of values of K within the vector of distances \code{s}
 for each point}
\item{s}{\code{s}}
}
\details{
The K function is defined as the expected number of further points
within a distance s of an arbitrary point, divided by the
overall density of the points. 
In practice an edge-correction is required to avoid biasing the
estimation due to non-recording of points outside the polygon.

The newstyle argument and khat object were introduced in collaboration 
with Thomas de Cornulier to permit the mapping of 
counts or khats for chosen distance values, as in 
\url{ftp://pbil.univ-lyon1.fr/pub/mac/ADE/ADE4/DocThemPDFUS/Thema81.pdf}, p.18.
}
\seealso{\code{\link{Kenv.csr}}}
\examples{
data(cardiff)
s <- seq(2,30,2)
plot(s, sqrt(khat(as.points(cardiff), cardiff$poly, s)/pi) - s,
 type="l", xlab="Splancs - polygon boundary", ylab="Estimated L",
 ylim=c(-1,1.5))
newstyle <- khat(as.points(cardiff), cardiff$poly, s, newstyle=TRUE)
str(newstyle)
newstyle
apply(newstyle$khats, 2, sum)
plot(newstyle)
}
\references{
Ripley, B.D. 1976 The second-order analysis of stationary point processes,
\emph{ J. Appl. Prob, 13} 255-266;
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{http://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

