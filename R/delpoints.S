# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

delpoints <- function(pts,add=FALSE)
{
	if(!add)pointmap(pts)
	cat("Use button 1 to delete points - button 2 when finished.\n")
	deld <- identify(pts[,1],pts[,2],atpen=FALSE,offset=0,labels=rep("X",length(pts[,1])))
	if(length(deld)!=0){
		pts <- pts[-deld,]
		}
	pts
}

