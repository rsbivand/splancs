# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

stsecal <- function(pts,times,poly,tlim,s,tm)
{
bigmat <- matrix(stvmat(pts,times,poly,tlim,s,tm),nrow=length(s)*length(tm))
diagon <- diag(bigmat)
semat <- sqrt(t(matrix(diagon,ncol=length(s))))
semat
}

