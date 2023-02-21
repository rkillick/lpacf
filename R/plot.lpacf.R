plot.lpacf <-
function (x, type = "line", lags = 1:min(as.integer(10 * 
    log10(nrow(x$lpacf))), ncol(x$lpacf) - 1), tcex = 1, lcol = 1, 
    llty = 1, the.time = NULL, plot.it = TRUE, xlab, ylab, ...) 
{
    if (missing(xlab)) 
        xlab <- "Time"
    if (min(lags) < 1)
	stop("Lag too small. Smallest lag is 1")
    nlags <- length(lags)
    ntime <- nrow(x$lpacf)
    if (max(lags) > ncol(x$lpacf)) 
        stop("Maximum lag is too high")
    if (length(lcol) == 1) 
        lcol <- rep(lcol, length(lags))
    if (length(llty) == 1) 
        llty <- rep(llty, length(lags))
    if (length(lcol) != length(lags)) 
        stop("Length of lcol vector has to be 1 or the same as the length of the lags vector")
    if (length(llty) != length(lags)) 
        stop("Length of llty vector has to be 1 or the same as the length of the lags vector")
    if (type == "line") {
        if (plot.it == TRUE) {
        	if (missing(ylab)) 
                  ylab <- "Partial Autocorrelation"
                plot(range(x$the.x), c(-1, 1), type = "n", 
                  xlab = xlab, ylab = ylab, ...)
                for (i in 1:nlags) {
                  lines(x$the.x, x$lpacf[, lags[i]], col = lcol[i], 
                    lty = llty[i])
                  pp <- seq(from = 1, to = ntime, length = 5)
                  text(x$the.x[pp], x$lpacf[pp, lags[i]], labels = lags[i], 
                    cex = tcex)
                }
            }
            ans <- x$lpacf[, lags]
            dimnames(ans) <- list(NULL, as.character(lags))
            return(invisible(ans))
    }
    else if (type == "persp") {
        m <- x$lpacf[, lags]
        zlab <- "Partial AC"
        if (plot.it == TRUE) {
            if (missing(ylab)) 
                ylab <- "Lag"
            persp(x = 1:ntime, y = lags, z = m[, lags], xlab = xlab, 
                ylab = ylab, zlab = zlab, ...)
        }
        ans <- m
        dimnames(ans) <- list(NULL, lags)
        return(ans)
    }
    else if (type == "acf") {
        if (is.null(the.time)) 
            stop("You have to specify a time point at which you wish to see the autocovariance/correlation. Specify the.time")
	time.mse <- (x$the.x - the.time)^2
	closix <- which(time.mse == min(time.mse))
	the.real.time <- x$the.x[closix]
        pacvals <- x$lpacf[closix, lags]
            if (missing(ylab)) 
                ylab <- "Partial Autocorrelation"
        if (plot.it == TRUE) {
            if (missing(xlab) || xlab== "Time") 
                xlab <- "Lag"
            plot(c(1, max(lags)), c(min(pacvals, 0), 1), type = "n", 
                xlab = xlab, ylab = ylab,
		sub=paste("At actual time: ", the.real.time), ...)
            segments(x0 = lags, y0 = 0, x1 = lags, y1 = pacvals)
            abline(h = 0)
        }
	names(pacvals) <- as.character(lags)
        return(invisible(pacvals))
    }
}
