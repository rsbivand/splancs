# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2001 by Roger S. Bivand

pip <- function(pts,poly,out=FALSE,bound=NULL, quiet=TRUE)
{
	inoutv <- inout(pts,poly,bound,quiet)
	if(!out){
		res <- pts[inoutv==TRUE,,drop=FALSE]
	}
	else{
		res <- pts[inoutv==FALSE,,drop=FALSE]
	}
	res
}

