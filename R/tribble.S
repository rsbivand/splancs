# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

tribble <- function(ccflag,vars=NULL,alphas=NULL,betas=NULL,rho
  ,which=1:length(alphas),covars=NULL,thetas=NULL,steps=NULL
  ,reqmin=0.001,icount=50,hessian=NULL)
{
#        library.dynam('splancs','tribble.o')
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
# scale distances
		disscale <- apply(vars,2,max)
		vars <- vars/matrix(disscale,ncol=ndvars,nrow=npt,byrow=TRUE)
	}
	if(is.null(covars) | is.null(thetas)){
		ncovars <- 0
	} else {
		covars <- as.matrix(covars)
		ncovars <- dim(covars)[2]
		if(ncovars != length(thetas))stop('Number of parameters inconsistent with covariate array')
		if(dim(covars)[1] != npt)stop('Number of covariates inconsistent with number of cases/controls')
# scale covariates
		covscale <- abs(apply(covars,2,max))
		covars <- covars/matrix(covscale,ncol=ncovars,nrow=npt,byrow=TRUE)
		
	}

#	allvars <- ifelse(ncovars != 0, cbind(vars,covars),
#             matrix(vars, ncol=ndvars, nrow=npt))
# fixing that ifelse thing gives an R tribble function that
# reproduces exactly the Splus output. (B.Rowlingson@lancaster.ac.uk
# 13 Feb 2002)
        if (ncovars != 0){
          allvars <- cbind(vars, covars)
        }else{
          allvars <- matrix(vars, ncol=ndvars, nrow=npt)
        }

	nallpars <- ndvars*2+ncovars + 1
	
        print(nallpars)
        storage.mode(allvars) <- "double"
	
	pstart <- c(alphas,betas,thetas,rho)
# initial steps - .5 for a,b,t, .1 for rho
	if(is.null(steps))steps <- c(rep(0.5,length=ndvars*2+ncovars),.1)

	l <- .Fortran('tribble',
		as.integer(ccflag),
		(allvars),
		as.integer(npt),
		as.integer(ndvars),
		as.integer(ndpars),
		as.integer(ncovars),
		as.integer(which),
		as.double(pstart),
		parfin=as.double(pstart),
		as.double(steps),
		as.double(reqmin),
		icode=as.integer(icount),
		kcode=as.integer(1),
		dlogl=as.double(1.00),
		PACKAGE="splancs")
	alphas <- betas <- thetas <- NULL
	if(ndpars!=0){
		alphas <- l$parfin[1:ndpars]
		betas <- (l$parfin[(ndpars+1):(ndpars*2)])/disscale
	}
	if(l$icode < 0)stop(paste("Error in minimisation algorithm code ",
          l$icode))
	rho <- l$parfin[length(l$parfin)]
	if(ncovars != 0){
	  thetas <- l$parfin[(1+ndpars*2):(ncovars+ndpars*2)]
	  thetas <- thetas/covscale
	}

	ncase <- sum(ccflag)
        ncont <- sum(1-ccflag)
        p <- ncase/(ncase+ncont)
        null.logl <- log(p)*ncase+log(1-p)*ncont

res <- list(alphas=alphas,betas=betas,thetas=thetas,rho=rho,logl=l$dlogl,
null.logl=null.logl,kcode=l$kcode,
call=match.call())
	class(res) <- "ribfit"

	if(!is.null(hessian)){
          warning("no hessian functions available")
#	  hess1 <- hessian.ribfit(res,ccflag=ccflag,vars=vars,covars=covars,which=which)
#	  hess2 <- hessian2.ribfit(res,ccflag=ccflag,vars=vars,covars=covars,which=which)
#	  res$hessian.1 <- hess1
#	  res$hessian.2 <- hess2
	}

	res
}
	


# Local Variables:
# mode:S
# S-temp-buffer-p:t
# End:
