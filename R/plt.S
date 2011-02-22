# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2002 by Roger S. Bivand,
# revised following suggestion by Brian Ripley

plt <- function(data, value)
{
	data <- sort(data)
	nd <- length(data)
	res <- numeric(length(value))
	for(iv in 1:length(value))
		res[iv] <- sum(as.integer(data < value[iv]))/nd
	res
}
