% Copyright Barry Rowlingson  and 
% Peter Diggle (c) 1991-3; https://www.maths.lancs.ac.uk/~rowlings/Splancs/
% R port: copyright 1998-2000 by Roger S. Bivand

\name{tribble}
\alias{tribble}
\title{
Diggle-Rowlingson Raised Incidence Model
}
\description{
Fits the Diggle-Rowlingson Raised Incidence Model.
}
\usage{
tribble(ccflag, vars=NULL, alphas=NULL, betas=NULL, rho, 
 which=1:length(alphas), covars=NULL, thetas=NULL, 
 steps=NULL, reqmin=0.001, icount=50, hessian=NULL)
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
Initial value of the \code{alpha} parameters.
}
\item{betas}{
Initial value of the \code{beta} parameters.
}
\item{rho}{
Initial value of the \code{rho} parameter.
}
\item{which}{
Defines the mapping from sources to parameters.
}
\item{covars}{
A matrix of covariates to be modelled as log-linear terms. The element 
\code{covars[i,j]} is the value of covariate \code{j} for case/control \code{i}.
}
\item{thetas}{
Initial values of covariate parameters.
}
\item{steps}{
Step sizes for the Nelder-Mead simplex algorithm.
}
\item{reqmin}{
Tolerance for simplex algorithm
}
\item{icount}{
Iteration count for simplex algorithm
}
\item{hessian}{by default NULL, any other value causes hessian to be computed and returned}
}
\value{
The return value is a list with many components, and class \code{ribfit}. 
\item{alphas}{A vector of the alpha parameters at the maximum}
\item{betas}{A vector of the beta values at the maximum}
\item{rho}{The value of rho at the maximum}
\item{logl}{The maximised log-likelihood}
\item{null.logl}{The null log-likelihood}
\item{call}{The function call to \code{tribble}}
For further information see Diggle and Rowlingson (1993).
}
\seealso{
\code{\link{triblik}}, \code{\link{dsquare}}
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

