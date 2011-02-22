# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

is.points <- function(p)
{
	is <- FALSE
	if(is.array(p))
		if(length(dim(p))==2)
			if(dim(p)[2]>=2) is <- TRUE
	is
}

