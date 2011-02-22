# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

print.ribfit <- function(x,...){
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	if(!is.null(x$kcode))cat(paste("Kcode = ",x$kcode,'\n',sep=''))

	if(!is.null(x$alpha)){
	  cat("\nDistance decay parameters:\n")
	  pmat <- cbind(x$alpha,x$beta)
	  dimnames(pmat) <- list(NULL,c('Alpha','Beta'))
	  print(pmat,...)
	}
	if(!is.null(x$theta)){
	  cat("\nCovariate parameters:\n")
	  print(x$theta,...)
	}
	cat(paste("\nrho parameter : ",x$rho,'\n\n',sep=''))

	cat(paste("     log-likelihood : ",x$logl,'\n',sep=''))
	cat(paste("null log-likelihood : ",x$null,'\n',sep=''))
	cat(paste("\n"))
	cat(paste("        D = 2(L-Lo) : ",2*(x$logl-x$null),'\n',sep=''))


}

# Local Variables:
# mode:S
# S-temp-buffer-p:t
# End:
