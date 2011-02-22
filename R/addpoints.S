# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

addpoints <- function(pts,plot=FALSE,quiet=FALSE)
{
	if(!plot & !missing(pts)){
		pointmap(pts)
	}

	if(missing(pts))pts <- NULL

	if(!quiet)cat("Use button 1 to enter new points - button 2 when finished.\n")
	newpts <- locator(type="p")
	pts <- rbind(pts,cbind(newpts$x,newpts$y))
	pts
}
