\name{summary.lpacf}
\alias{summary.lpacf}
\title{Print out summary information about a \code{lpacf} object.
}
\description{Print out summary information about a \code{lpacf} object.
}
\usage{
\method{summary}{lpacf}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The \code{lpacf} object you want to print out summary info for.
}
  \item{\dots}{Other arguments
}
}
\details{Prints out number of times that we have localized partial
	autocorrelation for and the number of lags computed. The localized
	partial autocorrelation is computed at a number of time points that
	might not be identical to the times in the original series, this
	function prints out the min and max of the range of times and
	an indicator of whether the whole series' localized pacf was computed.
	The bandwidth associated with spectral smoothing is printed and
	a note made of whether it was computed automatically or supplied
	as an earlier argument.
}
\value{None
}
\references{
Killick, R., Knight, M. I., Nason, G.P. and Eckley, I. A. (2020) The local partial autocorrelation function and some applications. \emph{Electron. J. Statist.} \bold{14} (2), 3268-3314. DOI: \href{https://doi.org/}{10.1214/20-EJS1748}.
}
\author{
Guy Nason
}

\seealso{\code{\link{lpacf}}, \code{\link{print.lacf}}
}
\examples{
#
# See example for print.lacf
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
