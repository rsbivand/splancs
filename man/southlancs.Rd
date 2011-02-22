\name{southlancs}
\title{Cancer cases in Chorley-Ribble}
\usage{data(southlancs)}
\alias{southlancs}
\alias{southlancs.pts}
\alias{southlancs.cc}
\alias{old.incinerator}
\alias{southlancs.bdy}
\format{A data frame with 974 observations
  \tabular{rlll}{
    [,1]  \tab x  \tab numeric  \tab grid eastings (metres)\cr
    [,2]  \tab y   \tab numeric  \tab grid northings (metres)\cr
    [,3]  \tab cc  \tab numeric  \tab case/control, lung=0, larynx=1\cr
  }
}
\source{Diggle, Gatrell and Lovett, 1990, - Bailey and Gatrell 1995, ch. 3.
}
\description{Locations of cases of cancer of lung and larynx in
Chorley-Ribble, Lancashire. The data set is split into a points object \code{southlancs.pts} and a
case/control 0/1 vector \code{southlancs.cc}. There are 917 controls and 57
cases in this data set - these numbers differ from 978 and 58 in Diggle (1990)
and Diggle and Rowlingson (1994). The data set also includes the approximate
location of an old incinerator \code{old.incinerator}, as well as
\code{southlancs.bdy}, the study area boundary.

}
\references{Bailey and Gatrell 1995, ch. 3; Diggle, P. (1990) A point process
modelling approach to raised incidence of a rare phenomenon in the viscinity
of a prespecified point. Journal of the Royal Statistical Society, A, 153,
349-362; Diggle, P. and Rowlingson, B. (1994) A conditional approach to point
process modelling of elevated risk. Journal of the Royal Statistical Society,
A, 157, 433-440.

}
\examples{
data(southlancs)
op <- par(mfrow=c(2,1))
pointmap(southlancs.pts[southlancs.cc == 0,])
pointmap(old.incinerator, add=TRUE, col="red", pch=19)
title("Lung cancer controls")
pointmap(southlancs.pts[southlancs.cc == 1,])
pointmap(old.incinerator, add=TRUE, col="red", pch=19)
title("Larynx cancer cases")
par(op)
polymap(southlancs.bdy,border="grey")
contour(kernel2d(southlancs.pts[southlancs.cc == 0,], 
	southlancs.bdy, h=500, nx=100, ny=100), nlevels=20, 
	add=TRUE,drawlabels=FALSE)
pointmap(southlancs.pts[southlancs.cc == 1,], add=TRUE, pch=19,
	 col="green")
pointmap(old.incinerator, add=TRUE, pch=19, col="red")
title(xlab="h=500, quartic kernel")
title("Density map of control, green case points, red old incinerator")

}
\keyword{datasets}
