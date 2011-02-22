# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2009 by Roger S. Bivand


stmctest <- function(pts,times,poly,tlimits,s,tt,nsim,quiet=FALSE,
   returnSims=FALSE)
{
	stats <- NULL
	ntpts <- npts(pts)
	kstsim <- stkhat(pts,times,poly,tlimits,s,tt)
	resids <- kstsim$kst-outer(kstsim$ks,kstsim$kt)
        t0 <- sum(resids)
        if (returnSims) {
            obs <- kstsim
            sims <- vector(mode="list", length=nsim)
        }
	for(isim in (1:nsim))
	{
		if(!quiet)cat("Doing simulation ", isim, "\n")
		times <- sample(times)
		kstsim <- stkhat(pts,times,poly,tlimits,s,tt)
		resids <- kstsim$kst-outer(kstsim$ks,kstsim$kt)
		stats <- c(stats,sum(resids))
		rdat <- rank(c(t0, stats))[1]
		if(!quiet)cat("Data ranks ",rdat," of ",isim,"\n")
                if (returnSims) sims[[isim]] <- kstsim
	}
	res <- list(t0=t0, t=stats)
        if (returnSims) {
            attr(res, "obs") <- obs
            attr(res, "sims") <- sims
        }
        res
}


# Local Variables:
# mode:S
# S-temp-buffer-p:t
# End:
