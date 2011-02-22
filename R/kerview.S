# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

"kerview"<-
function(pts,times, k3, map = TRUE, addimg=TRUE, ncol=1)
{
        pts3 <- cbind(pts,times)
        zgr <- k3$zgr
	old.par <- par()
	im.dev <- dev.cur()
	par(mfrow=c(ncol,ncol))
	if(ncol!=1)addimg <- FALSE
	brks <- quantile(k3$v, seq(0,1,0.05))
	cols <- heat.colors(length(brks)-1)
	image(k3$xgr, k3$ygr, k3$v[,  , 1], breaks=brks, col=cols, asp=1)
	map.dev <- -999
	if(map) {
		pbox <- bboxx(bbox(pts3))
		dstring <- paste(names(dev.cur()), 
			"(\"-geometry 500x500-5+5\")", sep = "")
		eval(parse(text = dstring))	
	#  dev.copy(eval(parse(text=names(dev.cur()))),'-geometry=500x500-5+5')
		map.dev <- dev.cur()
		on.exit(if(map.dev != -999) dev.off(map.dev))
		par(pty = "s")
		par(mfrow=c(ncol,ncol))
#		par(mai=c(0,0,0,0))
		pointmap(pts3, asp=1)
		tmin <- min(pts3[, 3])
		tmax <- max(pts3[, 3])
		tbins <- (zgr[1:(length(zgr) - 1)] + zgr[2:length(zgr)])/2
	#  if(tbins[1] < tmin)tbins <- c(tmin,tbins)
#  if(tbins[length(tbins)] > tmax)tbins <- c(tbins,tmax)
		pbins <- cut(pts3[, 3], tbins)
	}
#dev.copy(eval(parse(text=names(dev.cur()))),'-geometry=500x250')
	dstring <- paste(names(dev.cur()), "(\"-geometry 700x300+0-5\")", sep
		 = "")
	eval(parse(text = dstring))
	hi.dev <- dev.cur()
	on.exit({
		dev.off(hi.dev)
		if(map.dev != -999)
			dev.off(map.dev)
	}
	)
	hist(pts3[, 3])
	repeat {
		tsl <- locator(1)
		if(length(tsl) == 0)
			break
		tsl <- tsl$x
		bdist <- abs(tsl - k3$zgr)
		bin <- (1:length(k3$zgr))[bdist == min(bdist)]
		bin <- bin[1]
		dev.set(im.dev)
		image(k3$xgr, k3$ygr, k3$v[,  , bin], add = addimg,
		      breaks=brks, col=cols, asp=1)
		if(map) {
			dev.set(map.dev)
#			ptsin <- pts3[pbins == bin,]
		       ptsin <- pts3[abs(pts3[,3]-k3$zgr[bin]) < k3$hz,]
			if(length(ptsin != 0)) {
				pointmap(pbox, type = "n", asp=1)
				pointmap(ptsin, add = TRUE)
				title(paste("time = ", format(k3$zgr[bin])))
			}
		}
		dev.set(hi.dev)
	}
#dev.off(hi.dev)
	dev.set(im.dev)
}


# Local Variables:
# mode:S
# S-temp-buffer-p:t
# End:
