% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{print.ribfit}
\alias{print.ribfit}
\title{
Display the fit from tribble()
}
\description{
Display the fit from \code{tribble}
}
\usage{
\method{print}{ribfit}(x, ...)
}
\arguments{
\item{x}{
An object returned from \code{tribble}
}
\item{...}{optional arguments to pass through to \code{print()}}
}
\details{
 The parameter estimates and log-likelihood for the raised incidence 
model are displayed. The likelihood ratio, D = 2*(L-Lo), is also given.
 This function is called whenever \code{print} operates on an object
with class \code{ribfit}.
}
\seealso{
\code{\link{tribble}}
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

