\name{plot.lpacf}
\alias{plot.lpacf}
\title{Plot localized partial autocorrelation information in an
	\code{lpacf} object. 
}
\description{Plots information contained within a \code{lpacf} object.
	Plot arrangement and options are similar to that in \code{plot.lacf}
	in the \code{locits} package. 
}
\usage{
\method{plot}{lpacf}(x, type = "line", lags = 1:min(as.integer(10 * log10(nrow(x$lpacf))),
    ncol(x$lpacf) - 1), tcex = 1, lcol = 1, llty = 1, the.time = NULL, plot.it = TRUE, 
    xlab, ylab, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The \code{lpacf} class object you wish to plot.
}
  \item{type}{The type of plot you want. This can be \code{"line"} where
	each partial autocorrelation is plotted as a line over time, or
	\code{"persp"} where the partial autocorrelation across time and lag
	is plotted or \code{"acf"} where a regular-style plot is drawn,
	like the one produced by the standard R \code{acf} function, and
	this is the partial autocorrelation fixed at a given point in time.
}
  \item{lags}{A vector of integers containing the lags you wish to show.
	Note that, unlike regular autocorrelation, the smallest lag is lag one.
}
  \item{tcex}{For the \code{type="line"} plot integers are displayed along
	the lines of their corresponding partial autocorrelations. This
	parameter controls the scaling of the lag numbers.
}
  \item{lcol}{A vector of colors, the same length as the \code{lags} vector
	which controls the colour of each line corresponding to each partial
	autocorrelation line.
}
  \item{llty}{As for \code{lcol} but the line types of each line.
}
  \item{the.time}{A time has to be specified for the \code{type="acf"} plot,
	as this argument supplies it.
}
  \item{plot.it}{If \code{TRUE} a plot is produced. If \code{FALSE} then
	no plot is produced but the function still executes and produces
	the same (invisible) output.
}
  \item{xlab}{A label for the x-axis.
}
  \item{ylab}{A label for the y-axis
}
  \item{\dots}{Other arguments to plot.
}
}
\details{Produces a graphical representation of localized
	partial autocorrelation.
}
\value{The localized partial autocorrelation values are returned. 
	Essentially the ones that are, or would have been, plotted are
	returned. The lags can be selected using the \code{lags} argument.
	A matrix is returned: each row corresponds to a time point, each
	column corresponds to a lag. All time points are returned.
	Only the lags specified in the \code{lags} component are returned.
	The \code{dimnames} component indicates which lags were returned.
}
\references{
Killick, R., Knight, M. I., Nason, G.P. and Eckley, I. A. (2020) The local partial autocorrelation function and some applications. \emph{Electron. J. Statist.} \bold{14} (2), 3268-3314. DOI: \href{https://doi.org/}{10.1214/20-EJS1748}.
}
\author{
Guy Nason
}
\note{
This function was adapted from \code{plot.lacf} from the \code{locits}
package.
}


\seealso{
\code{\link{lpacf}}
}
\examples{
#
# Generate a test series
#
x.test <- tvar2sim()
#
# Compute its localized partial autocorrelation
#
x.lpacf <- lpacf(x.test)
#
# Perform a line plot of the localized partial autocorrelation of x
#	draw the lags in colours 1 thru 5
#
#
plot(x.lpacf, lags=1:5, lcol=1:5)
#
# Now produce the same plot, but omit lag 3
#
plot(x.lpacf, lags=c(1,2,4,5), lcol=c(1,2,4,5))
#
# Now plot localized autocovariance around time 175 using the regular acf
# style plot. 
#
plot(x.lpacf, type="acf", the.time=175)
}
\keyword{hplot}
\keyword{ts}
