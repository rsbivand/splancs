# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

Kenv.label <- function(pts1,pts2,poly,nsim,s,quiet=FALSE)
{
	kmax <- rep(0,length=length(s))
	kmin <- rep(1.0E34,length=length(s))
	
	for(isim in (1:nsim)){
		if(!quiet)cat('Doing labelling ',isim,'/',nsim,'\n')
		labpts <- rLabel(pts1,pts2)
		k1sim <- khat(labpts[[1]],poly,s)
		k2sim <- khat(labpts[[2]],poly,s)
		diffk <- k1sim-k2sim
		kmax <- pmax(kmax,diffk)
		kmin <- pmin(kmin,diffk)
	}
	list(lower=kmin,upper=kmax)
}

