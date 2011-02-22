# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

dsquare <- 
function(pts, srcs, namepref="d")
{
  ptsxy <- as.points(pts)
  nameso <- names(pts)
	d <- NULL
	nsrcs <- npts(srcs)
	for(i in 1:nsrcs) {
		d <- cbind(d, (srcs[i, 1] - ptsxy[,1])^2 + (srcs[i, 2] - 
			ptsxy[,2])^2)
		nameso <- c(nameso,paste(namepref,i,sep=''))
	}
	res <- d
        colnames(res) <- nameso
        res
}

# Local Variables:
# mode:S
# S-temp-buffer-p:t
# End:
