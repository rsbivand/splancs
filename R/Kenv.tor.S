# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

Kenv.tor <- function(pts1,pts2,poly,nsim,s,quiet=FALSE)
{
	rect <- bboxx(bbox(poly))
	kmax <- rep(0, length = length(s))
	kmin <- rep(9.999999999999999e+33, length = length(s))
	for(isim in 1:nsim){
		if(!quiet)cat("Doing shift ", isim, "/", nsim, "\n")
		pts2 <- rtor.shift(pts2,rect)
		k12sim <- k12hat(pts1,pts2,poly,s)
		kmax <- pmax(kmax,k12sim)
		kmin <- pmin(kmin,k12sim)
	}
	list(lower=kmin,upper=kmax)
}

# Copyright (C) 2003 Stephen Eglen <stephen@inf.ed.ac.uk>

Kenv.tor1 <- function (pts1, pts2, poly, nsim, s, quiet = FALSE) 
{
  ## This function adapted from Kenv.tor in splancs library.
  ## This allows us to compute u for each simulation.
  
  rect <- bboxx(bbox(poly))
  kmax <- rep(0, length = length(s))
  kmin <- rep(1e+34, length = length(s))

  theor <- pi * (s^2)

  ## The equation for computing u, score.k, could be passed
  ## as an argument to this function.
  
  ## Equation for u from Diggle (1986)
  ##score.k <- function(k) { sum( ( sqrt(k) - sqrt(theor))^2) }

  ## Equation for u from Diggle (p59, 2002) 
  score.k <- function(k) { sum( (1/(s^2))* ( k - theor)^2) }
  
  u <- rep(0, length=nsim+1)
  ## Calculate u for real K12 function.
  k12.real <-  k12hat(pts1, pts2, poly, s)
  u[1] <- score.k(k12.real)
  
  ## Keep a copy of all the K functions; first row is the K function of
  ## the real data set, then subsequent rows are the simulations.
  ## There is no need to keep this, although I find it useful for
  ## further processing.
  ksim <- matrix(0, nrow=nsim+1, ncol=length(k12.real))
  ksim[1,] <- k12.real
  
  for (isim in 1:nsim) {
    if (!quiet) 
      cat("Doing shift ", isim, "/", nsim, "\n")
    pts2 <- rtor.shift(pts2, rect)
    k12sim <- k12hat(pts1, pts2, poly, s)
    kmax <- pmax(kmax, k12sim)
    kmin <- pmin(kmin, k12sim)
    ksim[isim+1,] <- k12sim
    u[isim+1] <- score.k(k12sim)
  }
  
  list(lower = kmin, upper = kmax,
       real=k12.real, u=u,
       ksim=ksim,
       rank=rank(u)[1] )
}


