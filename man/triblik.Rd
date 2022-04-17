% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{triblik}
\alias{triblik}
\title{
Log-likelihood for the Diggle-Rowlingson raised incidence model
}
\description{
Calculates the log-likelihood for the Diggle-Rowlingson raised incidence model.
}
\usage{
triblik(ccflag, vars=NULL, alphas=NULL, betas=NULL, rho, 
 which=1:length(alphas), covars=NULL, thetas=NULL)
}
\arguments{
\item{ccflag}{
Case-control flag : a vector of ones and zeroes.
}
\item{vars}{
A matrix where \code{vars[i,j]} is the distance squared from point \code{i} 
to source \code{j}.
}
\item{alphas}{
The \code{alpha} parameters.
}
\item{betas}{
The \code{beta} parameters.
}
\item{rho}{
The \code{rho} parameter.
}
\item{which}{
Defines the mapping from sources to parameters.
}
\item{covars}{
A matrix of covariates to be modelled as log-linear terms. The element 
\code{covars[i,j]} is the value of covariate \code{j} for case/control \code{i}.
}
\item{thetas}{
The covariate parameters.
}}
\value{
The log-likelihood for the given parameters and the given distances and
optional covariates is returned.
}
\seealso{
\code{\link{tribble}}, \code{\link{dsquare}}
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

