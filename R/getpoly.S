# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

getpoly <- function(quiet=FALSE)
{
	if(!quiet){
		cat("Enter points with button 1\n")
		cat("Finish with button 2\n")
		cat("Don't try to join the polygon up - it is done for you.\n")
	}
	raw <- locator(type = "l")
	resloc <- cbind(raw$x, raw$y)
	lines(c(raw$x[1], raw$x[length(raw$x)]),
	      c(raw$y[1], raw$y[length(raw$y)]))
	resloc
}
