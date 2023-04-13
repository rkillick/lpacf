\name{lpacf.Epan}
\alias{lpacf.Epan}
%\alias{Epanechnikov}
\title{
Calculates the lpacf (local partial autocorrelation function) using an Epanechnikov kernel spectral smoother.
}
\description{
This calculates the local partial autocorrelation function for data \code{x} based on spectral smoothing using the Epanechnikov kernel.  Up to \code{lag.max} lags are calculated.
}
\usage{
lpacf.Epan(x, binwidth, lag.max=NULL, filter.number=10, family="DaubLeAsymm", 
    smooth.dev=var, AutoReflect=TRUE, tol=0.1, maxits=5, ABBverbose=0, lapplyfn=lapply,
    allpoints=FALSE, verbose=FALSE, ...)
}
\arguments{
  \item{x}{The time series you wish to analyze
}
  \item{binwidth}{The bandwidth for the spectral smoothing. If the argument is missing or zero then the bandwidth
  is chosen automatically using the \code{AutoBestBW.Epan} function. The smoothing is done using the Epanechnikov kernel in the \code{Epanechnikov} function.
}
  \item{lag.max}{The maximum lag of pacf required. If this argument is
	\code{NULL} then the maximum lag is computed internally to
	be \code{floor(10*log10(n))} where n is the length of the series.
	E.g. if n=512 then the \code{lag.max} will be set to 27.
}
  \item{filter.number}{The wavelet filter number for helping choose
	the smoothing \code{binwidth}, used only by \code{AutoBestBW}.
}
  \item{family}{The wavelet family for helping choose the smoothing
	\code{binwidth}, used only by \code{AutoBestBW}.
}
  \item{smooth.dev}{The type of deviance used by \code{AutoBestBW}
	to smooth an internal wavelet spectral estimate. The usual
	variance \code{var} is a good option, but something like a
	rescaled \code{mad} could be an alternative.
}
  \item{AutoReflect}{
	Mitigate periodic boundary conditions of wavelet transforms
	by reflecting time series about RHS end before taking 
	transforms (and is undone before returning the answer).	
	Setting this to be \code{TRUE} usually dramatically improves
	the results as the wavelet smoothing is performed using periodic
	transforms.
}
  \item{tol}{
	Tolerance for golden section search for the best bandwidth used
	only by \code{AutoBestBW}.
}
  \item{maxits}{
  Maximum number of iterations for the golden section search used by
	\code{AutoBestBW}.
}
  \item{ABBverbose}{
  If nonzero prints out informative messages about the progress
	of the golden section search. Higher integers produce more
	messages. Setting it to zero suppresses messages.
}
  \item{lapplyfn}{Function which applies a specified function to a vector.  lapply is default but this argument can be used to replace this by something more efficient, e.g. a parallel version such as \code{mclapply}. Note, if you
use \code{mclapply} then you should set the option \code{mc.cores} to something
sensible. E.g. I have a quad-core machine, so when I am using the machine
myself, alone, I set \code{options(mc.cores=4)}. If you have many cores on a
multicore machine, but with several users, you might want the number of
codes to be used by your code to be less than the max to be nice to others.

  }
  \item{allpoints}{The lpacf is calculated using a window centered on a time point, if allpoints=TRUE then the edges of the data are also estimated (where a smaller, non-centered binwidth is used). Note, if \code{allpoints=TRUE} you obtain
	an estimator for the whole length of the series, but it will be more
	variable nearer the ends.}
  \item{verbose}{If \code{TRUE} then some informative messages are printed,
	otherwise they're not}
  \item{\dots}{Other arguments for \code{\link[locits]{AutoBestBW}}.
}
}
\details{
Calculates the local partial autocorrelation function (lpacf) for 1,...,\code{lag.max} lags.   See paper in the references for more details.
NOTE: Often when local (windowed) estimates are created one assigns the estimated value to the central point in the window.  This is the approach we take here when calculating the lacv and lpacf.  This differs from the lpacf calculated for the \code{\link[forecastLSW]{forecast.lpacf}} function which assign the estimated value to the last point in the window.

The function works by using the regular R \code{pacf} applied to carefully
chosen windows of the original series. Note: code in the
\code{forecast.lpacf} from the \code{forecastLSW} package can compute a slightly different version of the
localized pacf using a wavelet method.
}
\value{An object of class \code{lpacf}. This is a list with the following
components:
\item{the.x}{the time coordinates of the localized partial autocorrelations.
	Call the length of this vector n.
	These coordinates can be non-integer (usually at halves) even
	if the original times are at the integers (which is assumed by
	this function). You can think of these indexing the value of
	the localized partial autocorrelation centred at these locations.}
\item{lpacf}{The localized partical autocorrelations.
Matrix of dimension n x \code{lag.max} containing the lpacf for each time point at lags 1,...,\code{lag.max}.}
\item{the.vacc}{This is \code{TRUE} if \code{allpoints=TRUE} and \code{FALSE}
	otherwise.}
\item{the.x1}{If \code{allpoints=TRUE} this vector contains the 
	x coordinate (time coordinate) values of the left-hand end of
	the series that are computed using increasingly reduced numbers of
	data points.}
\item{the.x2}{As for \code{the.x1} but for the right-hand end of the series.}
\item{vacc}{A vector of length two containing the interval of time points for
	which the localized autocovariance is computed using the largest
	number of points in the calculation --- ie the full binwidth.}
\item{binwidth}{The smoothing binwidth that was used.}
\item{AutoBinWidth}{This is \code{TRUE} if the binwidth was automatically
	selected, and \code{FALSE} if it was not.}
}
\references{
Killick, R., Knight, M. I., Nason, G.P. and Eckley, I. A. (2020) The local partial autocorrelation function and some applications. \emph{Electron. J. Statist.} \bold{14} (2), 3268-3314. DOI: \href{https://www.doi.org/}{10.1214/20-EJS1748}. 
}
\author{
Guy Nason and Rebecca Killick
}

\seealso{
\code{\link{lpacf}}, \code{\link{lpacf.plot}}, \code{\link[forecastLSW]{forecast.lpacf}},
\code{\link{plot.lacf}}, \code{\link{print.lpacf}}, \code{\link{summary.lpacf}}
}
\examples{
\donttest{
# first generate a time-varying process
set.seed(1)
x=tvar2sim()

x.lpacf <- lpacf(x)

#
# There are two functions to plot lpacf class objects
#
# One is via the generic function plot which uses plot.lpacf
# the other is a bespoke function lpact.plot. We'll look at
# the generic function first which behaves similarly to the equivalent
# function in the locits package: plot.lacf
#
#
# This plot shows all of the localized partial autocovariances up to lag 27
# which is the default calculated value for this length of time series.
plot(x.lpacf)

#
# In the previous plot, maybe there were too many lags. So, let's restrict to
# five lags and colour them differently.
#
plot(x.lpacf, lags=1:5, lcol=1:5)
#
# By default, the lpacf is not computed for the whole time series range.
# Let's do it for all the points now, and replot.
#
x.lpacf.all <- lpacf(x, allpoints=TRUE)
plot(x.lpacf.all, lags=1:5, lcol=1:5)
#
#
# Suppose we wanted to look at the localized partial autocorrelation at a
# particular time point, using the regular acf-like plot. We can do this by:
#
# We will choose the time point to examine the localized pacf at as 150.
#
plot(x.lpacf, type="acf", the.time=150)

# calculate the lpacf
ans<-lpacf.Epan(x,lag.max=10,filter.number=2,family="DaubExPhase")

# then maybe plot it by lag
lpacf.plot(ans,atLag=1:10,atTime=100)
}
}
\keyword{ts}
\keyword{wavelet}
