% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{stsecal}
\alias{stsecal}
\title{
Standard error for space-time clustering
}
\description{
Computes the standard error for space-time clustering.
}
\usage{
stsecal(pts, times, poly, tlim, s, tm)
}
\arguments{
\item{pts}{
A set of points, as defined in Splancs.
}
\item{times}{
A vector of times, the same length as the number of points in \code{pts}
}
\item{poly}{
A polygon enclosing the points
}
\item{tlim}{
A vector of length 2 specifying the upper and lower temporal domain.
}
\item{s}{
A vector of spatial distances for the analysis
}
\item{tm}{
A vector of times for the analysis
}}
\value{
A matrix of dimension \code{[length(s),length(t)]} is returned. Element 
\code{[i,j]} is the standard error at \code{s[i],t[j]}.
See Diggle Chetwynd Haggkvist and Morris (1995) for details.
}
\seealso{
\code{\link{stkhat}}, \code{\link{stsecal}}, \code{\link{stvmat}}, \code{\link{stdiagn}}
}

\references{
Diggle, P., Chetwynd, A., Haggkvist, R. and Morris, S. 1995 Second-order analysis of space-time clustering. Statistical Methods in Medical Research, 4, 124-136;Bailey, T. C. and Gatrell, A. C. 1995, Interactive spatial data analysis. Longman, Harlow, pp. 122-125;
Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis
code in S-Plus.  Computers and Geosciences, 19, 627-655;
the original sources can be accessed at:
\url{https://www.maths.lancs.ac.uk/~rowlings/Splancs/}. See also Bivand, R. and
Gebhardt, A. 2000 Implementing functions for spatial statistical analysis
using the R language. Journal of Geographical Systems, 2, 307-317.
}
\examples{
example(stkhat)
bur1se <- stsecal(burpts, burkitt$t, burbdy, c(400, 5800),
 seq(1,40,2), seq(100, 1500, 100))
}
% Converted by Sd2Rd version 0.3-3.
\keyword{spatial}

