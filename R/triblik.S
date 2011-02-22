# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

triblik <- function(ccflag,vars=NULL,alphas=NULL,betas=NULL
,rho,which=1:length(alphas),covars=NULL,thetas=NULL)
{

	vars <- cbind(vars)  # make sure its a matrix
	nalphas <- length(alphas)
	nbetas <- length(betas)
	if(nalphas!=nbetas){
		stop("alphas and betas must be same length")
	}
	else
	{
		ndpars <- nalphas
	}
	npt <- length(ccflag)
	if(is.null(vars)){
		ndvars <- 0
		if(is.null(covars))stop("No variates found!!")
	}
	else
	{
		ndvars <- dim(vars)[2]
		if(dim(vars)[1] != npt)stop('Number of variates inconsistent with number of cases/controls')
		if(ndvars < ndpars)stop('More source parameters than variables!')
		if(length(which) != ndvars) stop(' The "which" parameter is not the same length as the number of source parameters')
		wsu <- sort(unique(which))
		if(length(wsu) != ndpars)stop('Not enough distinct  values in "which"')
		if(!any(wsu == 1:ndpars))stop('Invalid values in "which"')
	}
	if(is.null(covars) | is.null(thetas)){
		ncovars <- 0
	} else {
		covars <- as.matrix(covars)
		ncovars <- dim(covars)[2]
		if(ncovars != length(thetas))stop('Number of parameters inconsistent with covariate array')
		if(dim(covars)[1] != npt)stop('Number of covariates inconsistent with number of cases/controls')
	}

	allvars <- cbind(vars,covars)
	nallpars <- ndvars*2+ncovars + 1
	
        storage.mode(allvars) <- "double"
	
	pars <- c(alphas,betas,thetas,rho)

	l <- .Fortran('trblik',
		as.integer(ccflag),
		(allvars),
		as.integer(npt),
		as.integer(nallpars),
		as.integer(ndvars),
		as.integer(which),
		as.integer(ndpars),
		as.double(pars),
		dlogl=as.double(1.00),
		PACKAGE="splancs")
	l$dlogl
}
	


# Local Variables:
# mode:S
# S-temp-buffer-p:t
# End:
