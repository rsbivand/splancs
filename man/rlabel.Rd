% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{rLabel}
\alias{rLabel}
\title{
Randomly label two or more point sets
}
\description{
Randomly label two or more point sets. (function name changed from rlabel to rLabel to avoid collision with spatstat)
}
\usage{
rLabel(\dots)
}
\arguments{
\item{...}{
Any number of points data sets
}}
\value{
A list of points data sets. There are as many elements in the list as 
arguments.


}
\details{
The output data sets are a random labelling of the input data sets, i.e. 
all the points in the input data sets are randomly assigned to the output 
sets. The number of points in each output set is the same as its corresponding
input set.
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

