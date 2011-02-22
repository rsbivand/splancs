% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{k12hat}
\alias{k12hat}
\title{
Bivariate K-function
}
\description{
Calculates an estimate of the bivariate K-function
}
\usage{
k12hat(pts1,pts2,poly,s)
}
\arguments{
\item{pts1,pts2}{
Two points data sets
}
\item{poly}{
A polygon containing the points
}
\item{s}{
A vector of distances at which to estimate the K12 function
}}
\value{
A vector like \code{s} containing the value of K12hat at the points in \code{s}.


}
\details{
The bivariate K function is defined as the expected number of points of
pattern 1
within a distance s of an arbitrary point of pattern 2, divided by the
overall density of the points in pattern 1.  To estimate this function,
the approximately unbiased estimator given by Lotwick and Silverman
(1982) is used.
}
\examples{
data(okwhite)
data(okblack)
okpoly <- list(x=c(okwhite$x, okblack$x), y=c(okwhite$y, okblack$y))
plot(seq(5,80,5), sqrt(k12hat(as.points(okwhite), as.points(okblack), 
bboxx(bbox(as.points(okpoly))), seq(5,80,5))/pi) - seq(5,80,5), xlab="distance", 
ylab=expression(hat(L)[12]), ylim=c(-20,20), type="l")
}
\references{
Lotwick, H.W. and Silverman B.W. (1982) Methods for Analysing Spatial
Processes of Several types of Points.
\emph{J. R. Statist Soc} B44 406-13;
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{http://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}

\keyword{spatial}

