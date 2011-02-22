# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

rLabel <- function(...){
	arglist <- list(...)
	nptsets <- length(arglist)
	ntotpts <- 0
	allpts <- NULL
	newperm <- NULL
	for(ia in 1:nptsets){
		nptsia <- npts(arglist[[ia]])
		ntotpts <- ntotpts+nptsia
		allpts <- rbind(allpts,arglist[[ia]])
		newperm <- c(newperm,rep(ia,times=nptsia))
	}
	newperm <- sample(newperm)
	outlist <- list()
	for(ia in 1:nptsets){
		outlist[[ia]] <- allpts[newperm==ia,]
	}
	outlist
}

